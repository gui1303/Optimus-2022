% LAC Test IEA15MW_03:  IEA 15 MW + Realistic wind preview
% Purpose:
% Here, we use a realistic wind preview to compare the behaviour with the
% platform damper on and off
% Result:
% 
% Authors:
% David Schlipf, Feng Guo
% Copyright (c) 2022 Flensburg University of Applied Sciences, WETI

%% Setup
clearvars;
close all;
clc;
addpath('..\MatlabFunctions')
addpath('..\MatlabFunctions\AnalyticlModel')

% Seeds (can be adjusted, but will provide different results)
nSample             = 6;                        % [-]           number of stochastic turbulence field samples
Seed_vec            = [1:nSample];              % [-]           vector of seeds

% Parameters postprocessing (can be adjusted, but will provide different results)
t_start             = 10;                       % [-]           ignore data before for STD and spectra
nDataPerBlock       = 570*40;                     % [-]           data per block, here 2^14/80 s = 204.8 s, so we have a frequency resolution of 1/204.8 Hz = 0.0049 Hz  
vWindow             = hamming(nDataPerBlock);   % [-]           window for estimation
nFFT                = [];                       % [-]           number of FFT, default: nextpow2(nDataPerBlock); 
nOverlap            = [];                       % [-]           samples of overlap, default: 50% overlap

% Files (should not be be changed)
TurbSimExeFile      = 'TurbSim_x64.exe';
FASTexeFile         = 'openfast_x64.exe';
FASTmapFile         = 'MAP_x64.lib';
SimulationName      = 'IEA-15-240-RWT-UMaineSemiNoDamper';
TurbSimTemplateFile = 'TurbSim2aInputFileTemplateIEA15MW.inp';
if ~exist('TurbulentWind','dir')
    mkdir TurbulentWind
end
if ~exist('SimulationResults','dir')
    mkdir SimulationResultsWithDamper
end
if ~exist('SimulationResults','dir')
    mkdir SimulationResultsNoDamper
end
%% Preprocessing: generate turbulent wind field
    
% Copy the adequate TurbSim version to the example folder 
copyfile(['..\TurbSim\',TurbSimExeFile],['TurbulentWind\',TurbSimExeFile])
    
% Generate all wind fields
for iSample = 1:nSample        
    Seed                = Seed_vec(iSample);
    TurbSimInputFile  	= ['TurbulentWind\URef_18_Seed_',num2str(Seed,'%02d'),'.ipt'];
    TurbSimResultFile  	= ['TurbulentWind\URef_18_Seed_',num2str(Seed,'%02d'),'.wnd'];
    if ~exist(TurbSimResultFile,'file')
        copyfile([TurbSimTemplateFile],TurbSimInputFile)
        ManipulateTXTFile(TurbSimInputFile,'MyRandSeed1',num2str(Seed));% adjust seed
        dos(['TurbulentWind\',TurbSimExeFile,' ',TurbSimInputFile]);
    end
end
    
% Clean up
delete(['TurbulentWind\',TurbSimExeFile])

%% Processing: run simulations with platform damper

ManipulateTXTFile('ROSCO_v2d6.IN','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');


% Copy the adequate OpenFAST version to the example folder
copyfile(['..\OpenFAST modified controller\',FASTexeFile],FASTexeFile)
copyfile(['..\OpenFAST modified controller\',FASTmapFile],FASTmapFile)

% Simulate with all wind fields
for iSample = 1:nSample
    
    % Adjust the InflowWind file
    Seed                = Seed_vec(iSample);
    WindFileRoot        = ['TurbulentWind\URef_18_Seed_',num2str(Seed,'%02d')];
    ManipulateTXTFile('IEA-15-240-RWT-UMaineSemi_InflowFile.dat','MyFilenameRoot',WindFileRoot);
    
    % Run FB    
    FASTresultFile      = ['SimulationResultsWithDamper\URef_18_Seed_WithDamper',num2str(Seed,'%02d'),'.outb'];
    if ~exist(FASTresultFile,'file')    
        dos([FASTexeFile,' ',SimulationName,'.fst']);
        movefile([SimulationName,'.outb'],FASTresultFile)
    end
   
    % Reset the InflowWind file again
    ManipulateTXTFile('IEA-15-240-RWT-UMaineSemi_InflowFile.dat',WindFileRoot,'MyFilenameRoot');
end

% Clean up
delete(FASTexeFile)
delete(FASTmapFile)


%% Processing: run simulations without platform damper

ManipulateTXTFile('ROSCO_v2d6.IN','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');

% Copy the adequate OpenFAST version to the example folder
copyfile(['..\OpenFAST modified controller\',FASTexeFile],FASTexeFile)
copyfile(['..\OpenFAST modified controller\',FASTmapFile],FASTmapFile)

% Simulate with all wind fields
for iSample = 1:nSample
    
    % Adjust the InflowWind file
    Seed                = Seed_vec(iSample);
    WindFileRoot        = ['TurbulentWind\URef_18_Seed_',num2str(Seed,'%02d')];
    ManipulateTXTFile('IEA-15-240-RWT-UMaineSemi_InflowFile.dat','MyFilenameRoot',WindFileRoot);
    
    % Run FB    
    FASTresultFile      = ['SimulationResultsNoDamper\URef_18_Seed_NoDamper',num2str(Seed,'%02d'),'.outb'];
    if ~exist(FASTresultFile,'file')    
        dos([FASTexeFile,' ',SimulationName,'.fst']);
        movefile([SimulationName,'.outb'],FASTresultFile)
    end
   
    % Reset the InflowWind file again
    ManipulateTXTFile('IEA-15-240-RWT-UMaineSemi_InflowFile.dat',WindFileRoot,'MyFilenameRoot');
end

% Clean up
delete(FASTexeFile)
delete(FASTmapFile)

%% Postprocessing: evaluate data with damper
%S_PtfmPitch_FB_est = NaN(nSample);
for iSample = 1:nSample    

    % Load data
    Seed                = Seed_vec(iSample);
    
    FASTresultFile      = ['SimulationResultsWithDamper\URef_18_Seed_WithDamper',num2str(Seed,'%02d'),'.outb'];
    FB_damper                  = ReadFASTbinaryIntoStruct(FASTresultFile);

    % Plot time results
    %figure('Name',['Seed ',num2str(Seed)])
    %hold on; grid on; box on
    %plot(FB_damper.Time,       FB_damper.Wind1VelX);
    %ylabel('RotSpeed [rpm]');
    %xlabel('time [s]')

    % Estimate spectra
    Fs                                      = 80; % [Hz]  sampling frequenzy, same as in *.fst
   [S_BldPitch1_FBWithDamper_est(iSample,:),f_est]	= pwelch(detrend(FB_damper.BldPitch1  (FB_damper.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
   [S_RotSpeed_FBWithDamper_est(iSample,:),f_est]	= pwelch(detrend(FB_damper.RotSpeed  (FB_damper.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
   [S_PtfmPitch_FBWithDamper_est(iSample,:),f_est]	= pwelch(detrend(FB_damper.PtfmPitch  (FB_damper.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    
end


%% Postprocessing: evaluate data without damper
%S_PtfmPitch_FB_est = NaN(nSample);
for iSample = 1:nSample    

    % Load data
    Seed                = Seed_vec(iSample);
    
    FASTresultFile      = ['SimulationResultsNoDamper\URef_18_Seed_NoDamper',num2str(Seed,'%02d'),'.outb'];
    FB_no_damper        = ReadFASTbinaryIntoStruct(FASTresultFile);

    % Plot time results
    %figure('Name',['Seed ',num2str(Seed)])
    %hold on; grid on; box on
    %plot(FB_no_damper.Time,       FB_no_damper.Wind1VelX);
    %ylabel('RotSpeed [rpm]');
    %xlabel('time [s]')

    % Estimate spectra
    Fs                                      = 80; % [Hz]  sampling frequenzy, same as in *.fst
    [S_BldPitch1_FBNoDamper_est(iSample,:),f_est]	= pwelch(detrend(FB_no_damper.BldPitch1  (FB_no_damper.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    [S_RotSpeed_FBNoDamper_est(iSample,:),f_est]	= pwelch(detrend(FB_no_damper.RotSpeed  (FB_no_damper.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    [S_PtfmPitch_FBNoDamper_est(iSample,:),f_est]	= pwelch(detrend(FB_no_damper.PtfmPitch  (FB_no_damper.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    
end

%% Plot spectra
figure('Name','Simulation results - platform pitch')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_PtfmPitch_FBWithDamper_est,1),'-','Color',[0 0.4470 0.7410]);
p2 = plot(f_est ,mean(S_PtfmPitch_FBNoDamper_est,1),'-','Color',[0.8500 0.3250 0.0980]);

set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Platform Pitch [(rpm)^2Hz^{-1}]') 
legend('Platform damper ON','Platform damper OFF')
hold off 


figure('Name','Simulation results - rotor speed')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_RotSpeed_FBWithDamper_est,1),'-','Color',[0 0.4470 0.7410]);
p2 = plot(f_est ,mean(S_RotSpeed_FBNoDamper_est,1),'-','Color',[0.8500 0.3250 0.0980]);

set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Rotor Speed [(rpm)^2Hz^{-1}]') 
legend('Platform damper ON','Platform damper OFF')
hold off


figure('Name','Simulation results - blade pitch')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_BldPitch1_FBWithDamper_est,1),'-','Color',[0 0.4470 0.7410]);
p2 = plot(f_est ,mean(S_BldPitch1_FBNoDamper_est,1),'-','Color',[0.8500 0.3250 0.0980]);

set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Blade Pitch [(rpm)^2Hz^{-1}]') 
legend('Platform damper ON','Platform damper OFF')
hold off