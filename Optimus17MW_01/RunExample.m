% LAC Test IEA15MW_03:  IEA 15 MW + Realistic wind preview
% Purpose:
%
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

% Parameters postprocessing (can be adjusted, but will provide different results)
t_start             = 30;                       % [-]           ignore data before for STD and spectra
nDataPerBlock       = 600*40;                   % [-]           data per block, here 2^14/80 s = 204.8 s, so we have a frequency resolution of 1/204.8 Hz = 0.0049 Hz  
vWindow             = hamming(nDataPerBlock);   % [-]           window for estimation
nFFT                = [];                       % [-]           number of FFT, default: nextpow2(nDataPerBlock); 
nOverlap            = [];                       % [-]           samples of overlap, default: 50% overlap

vWindSpeed          = 18:1:18;
NumSim              = length(vWindSpeed);

% Files (should not be be changed)
FASTexeFile         = 'openfast_x64.exe';
FASTmapFile         = 'MAP_x64.lib';
SimulationName      = 'IEA-15-255-RWT-UMaineSemi';

if ~exist('SimulationResultsConstant','dir')
    mkdir SimulationResultsConstant
end


%% Processing: run simulations with gust

% Copy the adequate OpenFAST version to the example folder
copyfile(['..\OpenFAST modified controller\',FASTexeFile],FASTexeFile)
copyfile(['..\OpenFAST modified controller\',FASTmapFile],FASTmapFile)

% Run the simulations
for iSim = 1:NumSim
    % Change wind file
    WindSpeed         = vWindSpeed(iSim);
    ManipulateTXTFile('IEA-15-255-RWT_UMaineSemi_InflowFile.dat','SetWind                HWindSpeed',[num2str(vWindSpeed(iSim)),'                HWindSpeed']);
    FASTresultFile      = ['SimulationResultsConstant\URef_18_Constant_',num2str(WindSpeed),'.outb'];
    if ~exist(FASTresultFile,'file')    
        dos([FASTexeFile,' ',SimulationName,'.fst']);
        movefile([SimulationName,'.outb'],FASTresultFile)
    end     
    % read in data
    FB_Constant(iSim)    = ReadFASTbinaryIntoStruct(FASTresultFile);
    % Reset the InflowWind file again
    ManipulateTXTFile('IEA-15-255-RWT_UMaineSemi_InflowFile.dat',[num2str(vWindSpeed(iSim)),'                HWindSpeed'],'SetWind                HWindSpeed');
    
    % Plot time results
    figure('Name','Time results for wind speed')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).Wind1VelX,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Wind velocity [m/s]');
    xlabel('time [s]')
    xlim([0 630]);
    hold off

    figure('Name','Time results for rotor speed')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).RotSpeed,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Rotor speed [RPM]');
    xlabel('time [s]')
    xlim([0 630]);
    hold off

    figure('Name','Time results for platform pitch')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).PtfmPitch,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Platform pitch [deg]');
    xlabel('time [s]')
    xlim([0 630]);
    hold off

    figure('Name','Time results for blade pitch')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).BldPitch1,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Blade pitch [deg]');
    xlabel('time [s]')
    xlim([0 630]);
    hold off

    figure('Name','Time results for electric power')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).RotPwr/1000,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Electric power [MW]');
    xlabel('time [s]')
    xlim([0 630]);
    hold off
end 

delete(FASTexeFile)
delete(FASTmapFile)



    
    

 