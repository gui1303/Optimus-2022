% Authors:
% David Schlipf, Feng Guo
% Copyright (c) 2022 Flensburg University of Applied Sciences, WETI

%% Setup
clearvars;
close all;
clc;
addpath('..\MatlabFunctions')
addpath('..\MatlabFunctions\AnalyticlModel')

vWindSpeed          = [15];
NumSim              = length(vWindSpeed);
vGenPwr             = NaN(NumSim,1);

% Files (should not be be changed)
FASTexeFile         = 'openfast_x64.exe';
FASTmapFile         = 'MAP_x64.lib';
SimulationName      = 'IEA-15-255-RWT-UMaineSemi';

if ~exist('SimulationResultsConstant','dir')
    mkdir SimulationResultsConstant
end

%% Processing: run simulations with constant wind

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
    vGenPwr(iSim)        = mean(FB_Constant(iSim).GenPwr);
    % Reset the InflowWind file again
    ManipulateTXTFile('IEA-15-255-RWT_UMaineSemi_InflowFile.dat',[num2str(vWindSpeed(iSim)),'                HWindSpeed'],'SetWind                HWindSpeed');
%     
%     % Plot time results
%     figure('Name','Time results for wind speed')
%     title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
%     hold on; grid on; box on
%     plot(FB_Constant(iSim).Time,FB_Constant(iSim).Wind1VelX,'Color',[0.8500 0.3250 0.0980]);
%     ylabel('Wind velocity [m/s]');
%     xlabel('time [s]')
%     xlim([0 630]);
%     hold off

    figure('Name','Time results for rotor speed')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).RotSpeed,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Rotor speed [RPM]');
    xlabel('time [s]')
    ylim([4 12]);
    xlim([0 630]);
    hold off

    figure('Name','Time results for platform pitch')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).PtfmPitch,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Platform pitch [deg]');
    xlabel('time [s]')
    xlim([0 630]);
    ylim([0 12]);
    hold off

    figure('Name','Time results for blade pitch')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).BldPitch1,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Blade pitch [deg]');
    xlabel('time [s]')
    xlim([0 630]);
    ylim([5 25]);
    hold off

    figure('Name','Time results for electric power')
    title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
    hold on; grid on; box on
    plot(FB_Constant(iSim).Time,FB_Constant(iSim).GenPwr/1000,'Color',[0.8500 0.3250 0.0980]);
    ylabel('Electric power [MW]');
    xlabel('time [s]')
    xlim([0 630]);
    ylim([0 22]);
    hold off    
end 

delete(FASTexeFile)
delete(FASTmapFile)

%plot power curve
figure('Name','Power curve')
title(['Wind speed ', num2str(vWindSpeed(iSim)), ' m/s'])
plot(vWindSpeed,vGenPwr,'-o');
ylabel('Electric power [MW]');
xlabel('Wind speed [m/s]')





 