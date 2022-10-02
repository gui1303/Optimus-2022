% Copy of RunExample from the IEW15MW_02 from the 
% https://github.com/fengguoFUAS/Baseline-Lidar-assisted-Controller/commit/0d83d87b6c7b03013fc2ad6e12c4fe993419dbd2)
% Modified to a simpler version without lidar and FF control.

%% Setup
clearvars; 
close all; 
clc;
addpath('..\MatlabFunctions')

% Copy the adequate OpenFAST version to the example folder
FASTexeFile     = 'openfast_x64.exe';
FASTmapFile     = 'MAP_x64.dll';
SimulationName  = 'IEA-15-240-RWT-UMaineSemi';
copyfile(['..\OpenFAST\',FASTexeFile],FASTexeFile)
copyfile(['..\OpenFAST\',FASTmapFile],FASTmapFile)

%% Run FB
dos([FASTexeFile,' ',SimulationName,'.fst']);                       % run OpenFAST
movefile([SimulationName,'.outb'],[SimulationName,'_FB.outb'])      % store results

%% Clean up
delete(FASTexeFile)
delete(FASTmapFile)

%% Comparison
% read in data
FB              = ReadFASTbinaryIntoStruct([SimulationName,'_FB.outb']);

% Plot         
ScreenSize = get(0,'ScreenSize');
figure('Name','Simulation results','position',[.1 .1 .8 .8].*ScreenSize([3,4,3,4]))

MyAxes(1) = subplot(4,1,1);
hold on; grid on; box on
plot(FB.Time,       FB.Wind1VelX);
ylabel('[m/s]');

MyAxes(2) = subplot(4,1,2);
hold on; grid on; box on
plot(FB.Time,       FB.BldPitch1);
ylabel('BldPitch1 [deg]');

MyAxes(3) = subplot(4,1,3);
hold on; grid on; box on
plot(FB.Time,       FB.RotSpeed);
ylabel('RotSpeed [rpm]');

MyAxes(4) = subplot(4,1,4);
hold on; grid on; box on
plot(FB.Time,       FB.PtfmPitch);
ylabel('PtfmPitch [deg]');

xlabel('time [s]')
linkaxes(MyAxes,'x');