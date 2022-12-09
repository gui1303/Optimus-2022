% Authors:
% David Schlipf, Feng Guo
% Copyright (c) 2022 Flensburg University of Applied Sciences, WETI

%% Setup
clearvars;
close all;
clc;

% Files (should not be be changed)
FASTexeFile         = 'openfast_x64.exe';
FASTmapFile         = 'MAP_x64.lib';
SimulationName      = 'IEA-15-255-RWT-UMaineSemi';

if ~exist('SimulationResultsConstant','dir')
    mkdir SimulationResultsConstant
end


%% Processing: run simulations with constant wind

% Run the simulation
FASTresultFile      = 'SimulationResultsConstant\URef_18_Constant_.outb';
if ~exist(FASTresultFile,'file')    
    dos([FASTexeFile,' ',SimulationName,'.fst']);
    movefile([SimulationName,'.outb'],FASTresultFile)
end

% read in data
FB_Constant    = ReadFASTbinaryIntoStruct(FASTresultFile);
figure('Name','Time results for rotor speed')
hold on; grid on; box on
plot(FB_Constant.Time,FB_Constant.RotSpeed,'Color',[0.8500 0.3250 0.0980]);
ylabel('Rotor speed [RPM]');
xlabel('time [s]')
hold off

figure('Name','Time results for platform pitch')
hold on; grid on; box on
plot(FB_Constant.Time,FB_Constant.PtfmPitch,'Color',[0.8500 0.3250 0.0980]);
ylabel('Platform pitch [deg]');
xlabel('time [s]')
hold off

figure('Name','Time results for blade pitch')
hold on; grid on; box on
plot(FB_Constant.Time,FB_Constant.BldPitch1,'Color',[0.8500 0.3250 0.0980]);
ylabel('Blade pitch [deg]');
xlabel('time [s]')
hold off

figure('Name','Time results for electric power')
hold on; grid on; box on
plot(FB_Constant.Time,FB_Constant.RotPwr/1000,'Color',[0.8500 0.3250 0.0980]);
ylabel('Electric power [MW]');
xlabel('time [s]')
hold off


    
    

 