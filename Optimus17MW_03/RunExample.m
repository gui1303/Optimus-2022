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
t_start             = 30;                       % [-]           ignore data before for STD and spectra
nDataPerBlock       = 600*2;                   % [-]           data per block, here 2^14/80 s = 204.8 s, so we have a frequency resolution of 1/204.8 Hz = 0.0049 Hz  
vWindow             = hamming(nDataPerBlock);   % [-]           window for estimation
nFFT                = [];                       % [-]           number of FFT, default: nextpow2(nDataPerBlock); 
nOverlap            = [];                       % [-]           samples of overlap, default: 50% overlap

% Files (should not be be changed)
TurbSimExeFile      = 'TurbSim_x64.exe';
FASTexeFile         = 'openfast_x64.exe';
FASTmapFile         = 'MAP_x64.lib';
SimulationName      = 'IEA-15-255-RWT-UMaineSemi';
TurbSimTemplateFile = 'TurbSim2aInputFileTemplateIEA15MW.inp';
if ~exist('TurbulentWind','dir')
    mkdir TurbulentWind
end
if ~exist('SimulationResultsTurbulentWind','dir')
    mkdir SimulationResultsTurbulentWind
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

%% Processing: run simulations with turbulent wind

% Copy the adequate OpenFAST version to the example folder
copyfile(['..\OpenFAST modified controller\',FASTexeFile],FASTexeFile)
copyfile(['..\OpenFAST modified controller\',FASTmapFile],FASTmapFile)

% Simulate with all wind fields
for iSample = 1:nSample
    
    % Adjust the InflowWind file
    Seed                = Seed_vec(iSample);
    WindFileRoot        = ['TurbulentWind\URef_18_Seed_',num2str(Seed,'%02d')];
    ManipulateTXTFile('IEA-15-255-RWT_UMaineSemi_InflowFile.dat','MyFilenameRoot',WindFileRoot);
    
    % Run FB    
    FASTresultFile      = ['SimulationResultsTurbulentWind\URef_18_Seed_TurbulentWind',num2str(Seed,'%02d'),'.outb'];
    if ~exist(FASTresultFile,'file')    
        dos([FASTexeFile,' ',SimulationName,'.fst']);
        movefile([SimulationName,'.outb'],FASTresultFile)
    end
   
    % Reset the InflowWind file again
    ManipulateTXTFile('IEA-15-255-RWT_UMaineSemi_InflowFile.dat',WindFileRoot,'MyFilenameRoot');
end

% Clean up
delete(FASTexeFile)
delete(FASTmapFile)





%% Obtain number of sampling points of the simulation
Seed                = 1;   
FASTresultFile      = ['SimulationResultsTurbulentWind\URef_18_Seed_TurbulentWind',num2str(Seed,'%02d'),'.outb'];
FB_TurbulentWind    = ReadFASTbinaryIntoStruct(FASTresultFile);
SimSamplingPoints   = size(FB_TurbulentWind.Time,1);
clear Seed FASTresultFile FB_TurbulentWind

%% Postprocessing: evaluate data with turbulent wind
%S_PtfmPitch_FB_est = NaN(nSample);
%Initialize avg power vector
AvgGenPwr                = NaN(nSample,1);
vRotSpeedTurbulentWind   = NaN(SimSamplingPoints,nSample);
vBldPitchTurbulentWind   = NaN(SimSamplingPoints,nSample);
vPtfmPitchTurbulentWind  = NaN(SimSamplingPoints,nSample);
for iSample = 1:nSample    

    % Load data
    Seed                = Seed_vec(iSample);
    
    FASTresultFile      = ['SimulationResultsTurbulentWind\URef_18_Seed_TurbulentWind',num2str(Seed,'%02d'),'.outb'];
    FB_TurbulentWind    = ReadFASTbinaryIntoStruct(FASTresultFile);
    AvgGenPwr(iSample)  = mean(FB_TurbulentWind.GenPwr);

%     % Plot time results
    figure('Name','Time results for rotor speed')
    hold on; grid on; box on
    plot(FB_TurbulentWind.Time,FB_TurbulentWind.RotSpeed,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Rotor speed');
    xlabel('time [s]')
    xlim([0 630]);
    hold off
    
    figure('Name','Time results for platform pitch')
    hold on; grid on; box on
    plot(FB_TurbulentWind.Time,FB_TurbulentWind.PtfmPitch,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Platform pitch [deg]');
    xlabel('time [s]')
    xlim([0 630]);
    hold off
    
    figure('Name','Time results for blade pitch')
    hold on; grid on; box on
    plot(FB_TurbulentWind.Time,FB_TurbulentWind.BldPitch1,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Blade pitch [deg]');
    xlabel('time [s]')
    xlim([0 630]);
    hold off
    
    figure('Name','Time results for electric power')
    hold on; grid on; box on
    plot(FB_TurbulentWind.Time,FB_TurbulentWind.GenPwr,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Electric power [kW]');
    xlabel('time [s]')
    xlim([0 630]);
    hold off
    
    

    % Load data for all wind fields in one variable to make the average
    for iSamplePoints = 1:SimSamplingPoints
        vRotSpeedTurbulentWind(iSamplePoints,iSample)  = FB_TurbulentWind.RotSpeed(iSamplePoints);
        vBldPitchTurbulentWind(iSamplePoints,iSample)  = FB_TurbulentWind.BldPitch1(iSamplePoints);
        vPtfmPitchTurbulentWind(iSamplePoints,iSample) = FB_TurbulentWind.PtfmPitch(iSamplePoints);
    end
    % Estimate spectra
    Fs                                      = 40; % [Hz]  sampling frequenzy, same as in *.fst
   [S_BldPitch1_TurbulenWind(iSample,:),f_est]	= pwelch(detrend(FB_TurbulentWind.BldPitch1  (FB_TurbulentWind.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
   [S_RotSpeed_TurbulentWind(iSample,:),f_est]	= pwelch(detrend(FB_TurbulentWind.RotSpeed  (FB_TurbulentWind.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
   [S_PtfmPitch_TurbulentWind(iSample,:),f_est]	= pwelch(detrend(FB_TurbulentWind.PtfmPitch  (FB_TurbulentWind.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    
end


% %% Plot avg spectra
figure('Name','Simulation results - platform pitch')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_PtfmPitch_TurbulentWind,1),'-','Color',[0.8500 0.3250 0.0980]);
set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Platform Pitch [(deg)^2Hz^{-1}]') 
hold off 


figure('Name','Simulation results - rotor speed')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_RotSpeed_TurbulentWind,1),'-','Color',[0.8500 0.3250 0.0980]);
set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Rotor Speed [(rpm)^2Hz^{-1}]') 
hold off


figure('Name','Simulation results - blade pitch')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_BldPitch1_TurbulenWind,1),'-','Color',[0.8500 0.3250 0.0980]);
set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Blade Pitch [(deg)^2Hz^{-1}]') 
hold off

AvgPwr = mean(AvgGenPwr)