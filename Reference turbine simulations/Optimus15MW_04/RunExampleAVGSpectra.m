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
nDataPerBlock       = 600*40;                   % [-]           data per block, here 2^14/80 s = 204.8 s, so we have a frequency resolution of 1/204.8 Hz = 0.0049 Hz  
vWindow             = hamming(nDataPerBlock);   % [-]           window for estimation
nFFT                = [];                       % [-]           number of FFT, default: nextpow2(nDataPerBlock); 
nOverlap            = [];                       % [-]           samples of overlap, default: 50% overlap

% Files (should not be be changed)
TurbSimExeFile      = 'TurbSim_x64.exe';
FASTexeFile         = 'openfast_x64.exe';
FASTmapFile         = 'MAP_x64.lib';
SimulationName      = 'IEA-15-240-RWT-UMaineSemi';
TurbSimTemplateFile = 'TurbSim2aInputFileTemplateIEA15MW.inp';
if ~exist('TurbulentWind','dir')
    mkdir TurbulentWind
end
if ~exist('SimulationResultsOriginalGain','dir')
    mkdir SimulationResultsOriginalGain
end
if ~exist('SimulationResults0_5Gain','dir')
    mkdir SimulationResults0_5Gain
end
if ~exist('SimulationResults0Gain','dir')
    mkdir SimulationResults0Gain
end
if ~exist('SimulationResults1_5Gain','dir')
    mkdir SimulationResults1_5Gain
end
if ~exist('SimulationResults2Gain','dir')
    mkdir SimulationResults2Gain
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

%% Processing: run simulations with original gain of the PD

%ManipulateTXTFile('ROSCO_v2d6.IN','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');
ManipulateTXTFile('ROSCO_v2d6.IN','-18.69290         ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]','-9.34645000000       ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]'); 


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
    FASTresultFile      = ['SimulationResultsOriginalGain\URef_18_Seed_OriginalGain',num2str(Seed,'%02d'),'.outb'];
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



%% Processing: run simulations with 50% of the original gain

ManipulateTXTFile('ROSCO_v2d6.IN','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');
ManipulateTXTFile('ROSCO_v2d6.IN','-9.34645000000       ! Fl_Kp             - Nacelle pitching proportional feedback ','-4.673225000       ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]'); 


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
    FASTresultFile      = ['SimulationResults0_5Gain\URef_18_Seed_0_5Gain',num2str(Seed,'%02d'),'.outb'];
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


%% Processing: run simulations without platform damper ( 0 gain)

%ManipulateTXTFile('ROSCO_v2d6.IN','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');
ManipulateTXTFile('ROSCO_v2d6.IN','-4.673225000       ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]','0      ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]'); 
% Just in case, reset the gain to the original value


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
    FASTresultFile      = ['SimulationResults0Gain\URef_18_Seed_0Gain',num2str(Seed,'%02d'),'.outb'];
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


%% Processing: run simulations with 50% increased gain

%ManipulateTXTFile('ROSCO_v2d6.IN','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');
ManipulateTXTFile('ROSCO_v2d6.IN','0      ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]','-14.01967500      ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]'); 


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
    FASTresultFile      = ['SimulationResults1_5Gain\URef_18_Seed_1_5Gain',num2str(Seed,'%02d'),'.outb'];
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


%% Processing: run simulations with double gain

%ManipulateTXTFile('ROSCO_v2d6.IN','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');
ManipulateTXTFile('ROSCO_v2d6.IN','-14.01967500      ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]','-18.69290         ! Fl_Kp             - Nacelle pitching proportional feedback gain [s]');


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
    FASTresultFile      = ['SimulationResults2Gain\URef_18_Seed_2Gain',num2str(Seed,'%02d'),'.outb'];
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

%% Obtain number of sampling points of the simulation
Seed                = 1;   
FASTresultFile      = ['SimulationResultsOriginalGain\URef_18_Seed_OriginalGain',num2str(Seed,'%02d'),'.outb'];
FB_OriginalGain     = ReadFASTbinaryIntoStruct(FASTresultFile);
SimSamplingPoints   = size(FB_OriginalGain.Time,1);
clear Seed FASTresultFile FB_OriginalGain

%% Postprocessing: evaluate data with original gain
%S_PtfmPitch_FB_est = NaN(nSample);
%figure('Name','Time results - original gain')
vRotSpeed_OriginalGain   = NaN(SimSamplingPoints,nSample);
vBldPitch_OriginalGain   = NaN(SimSamplingPoints,nSample);
vPtfmPitch_OriginalGain  = NaN(SimSamplingPoints,nSample);
for iSample = 1:nSample    

    % Load data
    Seed                = Seed_vec(iSample);
    
    FASTresultFile      = ['SimulationResultsOriginalGain\URef_18_Seed_OriginalGain',num2str(Seed,'%02d'),'.outb'];
    FB_OriginalGain     = ReadFASTbinaryIntoStruct(FASTresultFile);

%     % Plot time results
%     hold on; grid on; box on
%     plot(FB_OriginalGain.Time,       FB_OriginalGain.RotSpeed);
%     ylabel('Rotor speed');
%     xlabel('time [s]')
%     xlim([0 200]);

    % Load data for all wind fields in one variable to make the average
    for iSamplePoints = 1:SimSamplingPoints
        vRotSpeed_OriginalGain(iSamplePoints,iSample)  = FB_OriginalGain.RotSpeed(iSamplePoints);
        vBldPitch_OriginalGain(iSamplePoints,iSample)  = FB_OriginalGain.BldPitch1(iSamplePoints);
        vPtfmPitch_OriginalGain(iSamplePoints,iSample) = FB_OriginalGain.PtfmPitch(iSamplePoints);
    end
    % Estimate spectra
    Fs                                      = 40; % [Hz]  sampling frequenzy, same as in *.fst
   [S_BldPitch1_OriginalGain_est(iSample,:),f_est]	= pwelch(detrend(FB_OriginalGain.BldPitch1  (FB_OriginalGain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
   [S_RotSpeed_OriginalGain_est(iSample,:),f_est]	= pwelch(detrend(FB_OriginalGain.RotSpeed  (FB_OriginalGain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
   [S_PtfmPitch_OriginalGain_est(iSample,:),f_est]	= pwelch(detrend(FB_OriginalGain.PtfmPitch  (FB_OriginalGain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    
end

%% Postprocessing: evaluate data with 0% of original gain
%S_PtfmPitch_FB_est = NaN(nSample);
vRotSpeed_0Gain   = NaN(SimSamplingPoints,nSample);
vBldPitch_0Gain   = NaN(SimSamplingPoints,nSample);
vPtfmPitch_0Gain  = NaN(SimSamplingPoints,nSample);
%figure('Name','Time results - zero gain')
for iSample = 1:nSample    

    % Load data
    Seed                = Seed_vec(iSample);
    
    FASTresultFile           = ['SimulationResults0Gain\URef_18_Seed_0Gain',num2str(Seed,'%02d'),'.outb'];
    FB_0Gain                 = ReadFASTbinaryIntoStruct(FASTresultFile);
% 
%     % Plot time results
%     hold on; grid on; box on
%     plot(FB_0Gain.Time,       FB_0Gain.RotSpeed);
%     ylabel('Rotor speed');
%     xlabel('time [s]')
%     xlim([0 200]);

    % Load data for all wind fields in one variable to make the average
    for iSamplePoints = 1:SimSamplingPoints
        vRotSpeed_0Gain(iSamplePoints,iSample)  = FB_0Gain.RotSpeed(iSamplePoints);
        vBldPitch_0Gain(iSamplePoints,iSample)  = FB_0Gain.BldPitch1(iSamplePoints);
        vPtfmPitch_0Gain(iSamplePoints,iSample) = FB_0Gain.PtfmPitch(iSamplePoints);
    end
    
    % Estimate spectra
    Fs                                      = 40; % [Hz]  sampling frequenzy, same as in *.fst
   [S_BldPitch1_0Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_0Gain.BldPitch1  (FB_0Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
   [S_RotSpeed_0Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_0Gain.RotSpeed  (FB_0Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
   [S_PtfmPitch_0Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_0Gain.PtfmPitch  (FB_0Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    
end


%% Postprocessing: evaluate data 50% of original gain
%S_PtfmPitch_FB_est = NaN(nSample);
vRotSpeed_0_5Gain   = NaN(SimSamplingPoints,nSample);
vBldPitch_0_5Gain   = NaN(SimSamplingPoints,nSample);
vPtfmPitch_0_5Gain  = NaN(SimSamplingPoints,nSample);
%figure('Name','Time results - 0.5 times the original gain')
for iSample = 1:nSample    

    % Load data
    Seed                = Seed_vec(iSample);
    
    FASTresultFile      = ['SimulationResults0_5Gain\URef_18_Seed_0_5Gain',num2str(Seed,'%02d'),'.outb'];
    FB_0_5Gain        = ReadFASTbinaryIntoStruct(FASTresultFile);

%     % Plot time results
%     %figure('Name',['Seed ',num2str(Seed)])
%     hold on; grid on; box on
%     plot(FB_0_5Gain.Time,       FB_0_5Gain.RotSpeed);
%     ylabel('Rotor speed');
%     xlabel('time [s]')
%     xlim([0 200])

    
    % Load data for all wind fields in one variable to make the average
    for iSamplePoints = 1:SimSamplingPoints
        vRotSpeed_0_5Gain(iSamplePoints,iSample)  = FB_0_5Gain.RotSpeed(iSamplePoints);
        vBldPitch_0_5Gain(iSamplePoints,iSample)  = FB_0_5Gain.BldPitch1(iSamplePoints);
        vPtfmPitch_0_5Gain(iSamplePoints,iSample) = FB_0_5Gain.PtfmPitch(iSamplePoints);
    end
    
    
    % Estimate spectra
    Fs                                      = 40; % [Hz]  sampling frequenzy, same as in *.fst
    [S_BldPitch1_0_5Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_0_5Gain.BldPitch1  (FB_0_5Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    [S_RotSpeed_0_5Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_0_5Gain.RotSpeed  (FB_0_5Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    [S_PtfmPitch_0_5Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_0_5Gain.PtfmPitch  (FB_0_5Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    
end

%% Postprocessing: evaluate data 1.5 times the original gain
%S_PtfmPitch_FB_est = NaN(nSample);
vRotSpeed_1_5Gain   = NaN(SimSamplingPoints,nSample);
vBldPitch_1_5Gain   = NaN(SimSamplingPoints,nSample);
vPtfmPitch_1_5Gain  = NaN(SimSamplingPoints,nSample);

%figure('Name','Time results - 1.5 times the original gain')
for iSample = 1:nSample    

    % Load data
    Seed                = Seed_vec(iSample);
    
    FASTresultFile      = ['SimulationResults1_5Gain\URef_18_Seed_1_5Gain',num2str(Seed,'%02d'),'.outb'];
    FB_1_5Gain        = ReadFASTbinaryIntoStruct(FASTresultFile);

%     % Plot time results
%     %figure('Name',['Seed ',num2str(Seed)])
%     hold on; grid on; box on
%     plot(FB_1_5Gain.Time,       FB_1_5Gain.RotSpeed);
%     ylabel('Rotor speed');
%     xlabel('time [s]')
%     xlim([0 200])


    % Load data for all wind fields in one variable to make the average
    for iSamplePoints = 1:SimSamplingPoints
        vRotSpeed_1_5Gain(iSamplePoints,iSample)  = FB_1_5Gain.RotSpeed(iSamplePoints);
        vBldPitch_1_5Gain(iSamplePoints,iSample)  = FB_1_5Gain.BldPitch1(iSamplePoints);
        vPtfmPitch_1_5Gain(iSamplePoints,iSample) = FB_1_5Gain.PtfmPitch(iSamplePoints);
    end
    

    % Estimate spectra
    Fs                                      = 40; % [Hz]  sampling frequenzy, same as in *.fst
    [S_BldPitch1_1_5Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_1_5Gain.BldPitch1  (FB_1_5Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    [S_RotSpeed_1_5Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_1_5Gain.RotSpeed  (FB_1_5Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    [S_PtfmPitch_1_5Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_1_5Gain.PtfmPitch  (FB_1_5Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    
end



%% Postprocessing: evaluate data 2 times the original gain
%S_PtfmPitch_FB_est = NaN(nSample);
vRotSpeed_2Gain   = NaN(SimSamplingPoints,nSample);
vBldPitch_2Gain   = NaN(SimSamplingPoints,nSample);
vPtfmPitch_2Gain  = NaN(SimSamplingPoints,nSample);

%figure('Name','Time results - 2 times the original gain')
for iSample = 1:nSample    

    % Load data
    Seed                = Seed_vec(iSample);
    
    FASTresultFile      = ['SimulationResults2Gain\URef_18_Seed_2Gain',num2str(Seed,'%02d'),'.outb'];
    FB_2Gain            = ReadFASTbinaryIntoStruct(FASTresultFile);

%     % Plot time results
%     %figure('Name',['Seed ',num2str(Seed)])
%     hold on; grid on; box on
%     plot(FB_2Gain.Time,       FB_2Gain.RotSpeed);
%     ylabel('Rotor speed');
%     xlabel('time [s]')
%     xlim([0 200])
    
    % Load data for all wind fields in one variable to make the average
    for iSamplePoints = 1:SimSamplingPoints
        vRotSpeed_2Gain(iSamplePoints,iSample)  = FB_2Gain.RotSpeed(iSamplePoints);
        vBldPitch_2Gain(iSamplePoints,iSample)  = FB_2Gain.BldPitch1(iSamplePoints);
        vPtfmPitch_2Gain(iSamplePoints,iSample) = FB_2Gain.PtfmPitch(iSamplePoints);
    end
    
        
    % Estimate spectra
    Fs                                      = 40; % [Hz]  sampling frequenzy, same as in *.fst
    [S_BldPitch1_2Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_2Gain.BldPitch1  (FB_2Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    [S_RotSpeed_2Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_2Gain.RotSpeed  (FB_2Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    [S_PtfmPitch_2Gain_est(iSample,:),f_est]	= pwelch(detrend(FB_2Gain.PtfmPitch  (FB_2Gain.Time  >t_start)),vWindow,nOverlap,nFFT,Fs);
    
end


%% Plot spectra
figure('Name','Simulation results - platform pitch')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_PtfmPitch_OriginalGain_est,1),'-','Color',[0 0.4470 0.7410]);
p2 = plot(f_est ,mean(S_PtfmPitch_0Gain_est,1),'-','Color',[0.4500 0.3250 0.0980]);
p3 = plot(f_est ,mean(S_PtfmPitch_0_5Gain_est,1),'-','Color',[0.8500 0.3250 0.0980]);
p4 = plot(f_est ,mean(S_PtfmPitch_1_5Gain_est,1),'-','Color',[0.3500 0.1250 0.0980]);
p5 = plot(f_est ,mean(S_PtfmPitch_2Gain_est,1),'-','Color',[0.6500 0.7250 0.0980]);

set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Platform Pitch [(deg)^2Hz^{-1}]') 
legend('Original gain','No platform damper','0.5 times the original gain','1.5 times the original gain','2 times the original gain')
hold off 


figure('Name','Simulation results - rotor speed')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_RotSpeed_OriginalGain_est,1),'-','Color',[0 0.4470 0.7410]);
p2 = plot(f_est ,mean(S_RotSpeed_0Gain_est,1),'-','Color',[0.4500 0.3250 0.0980]);
p3 = plot(f_est ,mean(S_RotSpeed_0_5Gain_est,1),'-','Color',[0.8500 0.3250 0.0980]);
p4 = plot(f_est ,mean(S_RotSpeed_1_5Gain_est,1),'-','Color',[0.3500 0.1250 0.0980]);
p5 = plot(f_est ,mean(S_RotSpeed_2Gain_est,1),'-','Color',[0.6500 0.7250 0.0980]);

set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Rotor Speed [(rpm)^2Hz^{-1}]') 
legend('Original gain','No platform damper','0.5 times the original gain','1.5 times the original gain','2 times the original gain')

hold off


figure('Name','Simulation results - blade pitch')

hold on; grid on; box on
p1 = plot(f_est ,mean(S_BldPitch1_OriginalGain_est,1),'-','Color',[0 0.4470 0.7410]);
p2 = plot(f_est ,mean(S_BldPitch1_0Gain_est,1),'-','Color',[0.4500 0.3250 0.0980]);
p3 = plot(f_est ,mean(S_BldPitch1_0_5Gain_est,1),'-','Color',[0.8500 0.3250 0.0980]);
p4 = plot(f_est ,mean(S_BldPitch1_1_5Gain_est,1),'-','Color',[0.3500 0.1250 0.0980]);
p5 = plot(f_est ,mean(S_BldPitch1_2Gain_est,1),'-','Color',[0.6500 0.7250 0.0980]);

set(gca,'Xscale','log')
set(gca,'Yscale','log')
xlabel('Frequency [Hz] ')
ylabel('Spectra Blade Pitch [(deg)^2Hz^{-1}]') 
legend('Original gain','No platform damper','0.5 times the original gain','1.5 times the original gain','2 times the original gain')
hold off


%% Average over all the six simulations
vBldPitchAvg  = NaN (nSample,1);
vBldPitchSTD  = NaN (nSample,1);
vRotSpeedAvg  = NaN (nSample,1);
vRotSPeedSTD  = NaN (nSample,1);
vPtfmPitchAvg = NaN (nSample,1);
vPtfmPitchSTD = NaN (nSample,1);
% 0 gain
for iSample = 1:nSample
    vBldPitchAvg(iSample)  = mean(vBldPitch_0Gain(:,iSample));
    vBldPitchSTD(iSample)  = std(vBldPitch_0Gain(:,iSample));
    vRotSpeedAvg(iSample)  = mean(vRotSpeed_0Gain(:,iSample));
    vRotSPeedSTD(iSample)  = std(vRotSpeed_0Gain(:,iSample));
    vPtfmPitchAvg(iSample) = mean(vPtfmPitch_0Gain(:,iSample));
    vPtfmPitchSTD(iSample) = std(vPtfmPitch_0Gain(:,iSample));
end
AvgValues.Gain0.BldPitchAvg   = mean(vBldPitchAvg);
AvgValues.Gain0.BldPitchSTD   = mean(vBldPitchSTD);
AvgValues.Gain0.RotSpeedAvg   = mean(vRotSpeedAvg);
AvgValues.Gain0.RotSpeedSTD   = mean(vRotSPeedSTD);
AvgValues.Gain0.PtfmPitchAvg  = mean(vPtfmPitchAvg);
AvgValues.Gain0.PtfmPitchSTD  = mean(vPtfmPitchSTD);

% 0.5 of original gain
for iSample = 1:nSample
    vBldPitchAvg(iSample)  = mean(vBldPitch_0_5Gain(:,iSample));
    vBldPitchSTD(iSample)  = std(vBldPitch_0_5Gain(:,iSample));
    vRotSpeedAvg(iSample)  = mean(vRotSpeed_0_5Gain(:,iSample));
    vRotSPeedSTD(iSample)  = std(vRotSpeed_0_5Gain(:,iSample));
    vPtfmPitchAvg(iSample) = mean(vPtfmPitch_0_5Gain(:,iSample));
    vPtfmPitchSTD(iSample) = std(vPtfmPitch_0_5Gain(:,iSample));
end
AvgValues.Gain0_5.BldPitchAvg   = mean(vBldPitchAvg);
AvgValues.Gain0_5.BldPitchSTD   = mean(vBldPitchSTD);
AvgValues.Gain0_5.RotSpeedAvg   = mean(vRotSpeedAvg);
AvgValues.Gain0_5.RotSpeedSTD   = mean(vRotSPeedSTD);
AvgValues.Gain0_5.PtfmPitchAvg  = mean(vPtfmPitchAvg);
AvgValues.Gain0_5.PtfmPitchSTD  = mean(vPtfmPitchSTD);

% Original gain
for iSample = 1:nSample
    vBldPitchAvg(iSample)  = mean(vBldPitch_OriginalGain(:,iSample));
    vBldPitchSTD(iSample)  = std(vBldPitch_OriginalGain(:,iSample));
    vRotSpeedAvg(iSample)  = mean(vRotSpeed_OriginalGain(:,iSample));
    vRotSPeedSTD(iSample)  = std(vRotSpeed_OriginalGain(:,iSample));
    vPtfmPitchAvg(iSample) = mean(vPtfmPitch_OriginalGain(:,iSample));
    vPtfmPitchSTD(iSample) = std(vPtfmPitch_OriginalGain(:,iSample));
end
AvgValues.GainOriginal.BldPitchAvg   = mean(vBldPitchAvg);
AvgValues.GainOriginal.BldPitchSTD   = mean(vBldPitchSTD);
AvgValues.GainOriginal.RotSpeedAvg   = mean(vRotSpeedAvg);
AvgValues.GainOriginal.RotSpeedSTD   = mean(vRotSPeedSTD);
AvgValues.GainOriginal.PtfmPitchAvg  = mean(vPtfmPitchAvg);
AvgValues.GainOriginal.PtfmPitchSTD  = mean(vPtfmPitchSTD);

% 1.5 the original gain
for iSample = 1:nSample
    vBldPitchAvg(iSample)  = mean(vBldPitch_1_5Gain(:,iSample));
    vBldPitchSTD(iSample)  = std(vBldPitch_1_5Gain(:,iSample));
    vRotSpeedAvg(iSample)  = mean(vRotSpeed_1_5Gain(:,iSample));
    vRotSPeedSTD(iSample)  = std(vRotSpeed_1_5Gain(:,iSample));
    vPtfmPitchAvg(iSample) = mean(vPtfmPitch_1_5Gain(:,iSample));
    vPtfmPitchSTD(iSample) = std(vPtfmPitch_1_5Gain(:,iSample));
end
AvgValues.Gain1_5.BldPitchAvg   = mean(vBldPitchAvg);
AvgValues.Gain1_5.BldPitchSTD   = mean(vBldPitchSTD);
AvgValues.Gain1_5.RotSpeedAvg   = mean(vRotSpeedAvg);
AvgValues.Gain1_5.RotSpeedSTD   = mean(vRotSPeedSTD);
AvgValues.Gain1_5.PtfmPitchAvg  = mean(vPtfmPitchAvg);
AvgValues.Gain1_5.PtfmPitchSTD  = mean(vPtfmPitchSTD);

% 2 time the original gain
for iSample = 1:nSample
    vBldPitchAvg(iSample)  = mean(vBldPitch_2Gain(:,iSample));
    vBldPitchSTD(iSample)  = std(vBldPitch_2Gain(:,iSample));
    vRotSpeedAvg(iSample)  = mean(vRotSpeed_2Gain(:,iSample));
    vRotSPeedSTD(iSample)  = std(vRotSpeed_2Gain(:,iSample));
    vPtfmPitchAvg(iSample) = mean(vPtfmPitch_2Gain(:,iSample));
    vPtfmPitchSTD(iSample) = std(vPtfmPitch_2Gain(:,iSample));
end
AvgValues.Gain2.BldPitchAvg  = mean(vBldPitchAvg);
AvgValues.Gain2.BldPitchSTD  = mean(vBldPitchSTD);
AvgValues.Gain2.RotSpeedAvg  = mean(vRotSpeedAvg);
AvgValues.Gain2.RotSpeedSTD  = mean(vRotSPeedSTD);
AvgValues.Gain2.PtfmPitchAvg = mean(vPtfmPitchAvg);
AvgValues.Gain2.PtfmPitchSTD = mean(vPtfmPitchSTD);

% Plot averages
avgSTDBladePitchVector = [AvgValues.Gain0.BldPitchSTD AvgValues.Gain0_5.BldPitchSTD AvgValues.GainOriginal.BldPitchSTD AvgValues.Gain1_5.BldPitchSTD AvgValues.Gain2.BldPitchSTD];
x = [0 0.5 1 1.5 2];
avgSTDPtfmPitchVector = [AvgValues.Gain0.PtfmPitchSTD AvgValues.Gain0_5.PtfmPitchSTD AvgValues.GainOriginal.PtfmPitchSTD AvgValues.Gain1_5.PtfmPitchSTD AvgValues.Gain2.PtfmPitchSTD];
avgSTDRotorSpeedVector = [AvgValues.Gain0.RotSpeedSTD AvgValues.Gain0_5.RotSpeedSTD AvgValues.GainOriginal.RotSpeedSTD AvgValues.Gain1_5.RotSpeedSTD AvgValues.Gain2.RotSpeedSTD];
figure 
plot(x,avgSTDBladePitchVector,'-o','Color',[0.8500 0.3250 0.0980])
title('Blade pitch standard deviation for different values of platform damper gain')
xlabel('Platform damper gain (multiple of the original ROSCO gain)')
ylabel('Standard deviation [deg]') 
figure
plot(x,avgSTDPtfmPitchVector,'-o','Color',[0.8500 0.3250 0.0980])
title('Platform pitch standard deviation for different values of platform damper gain')
xlabel('Platform damper gain (multiple of the original ROSCO gain)')
ylabel('Standard deviation [deg]') 
figure
plot(x,avgSTDRotorSpeedVector,'-o','Color',[0.8500 0.3250 0.0980])
title('Rotor speed standard deviation for different values of platform damper gain')
xlabel('Platform damper gain (multiple of the original ROSCO gain)')
ylabel('Standard deviation [RPM]') 