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
% if ~exist(FASTresultFile,'file')    
    dos([FASTexeFile,' ',SimulationName,'.fst']);
    movefile([SimulationName,'.outb'],FASTresultFile)
% end

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
plot(FB_Constant.Time,FB_Constant.GenPwr/1000,'Color',[0.8500 0.3250 0.0980]);
ylabel('Electric power [MW]');
xlabel('time [s]')
hold off


%% Plot         
ScreenSize = get(0,'ScreenSize');
figure('Name','Simulation results','position',[.1 .1 .8 .8].*ScreenSize([3,4,3,4]))
n = 5;

MyAxes(1) = subplot(n,1,1);
hold on; grid on; box on
plot(FB_Constant.Time,       FB_Constant.Wind1VelX);
ylabel('[m/s]');

MyAxes(2) = subplot(n,1,2);
hold on; grid on; box on
plot(FB_Constant.Time,       FB_Constant.BldPitch1, 'r');
ylabel('BldPitch1 [deg]');

MyAxes(3) = subplot(n,1,3);
hold on; grid on; box on
plot(FB_Constant.Time,       FB_Constant.GenTq, 'r');
ylabel('GenTq [kNm]');

MyAxes(4) = subplot(n,1,4);
hold on; grid on; box on
plot(FB_Constant.Time,       FB_Constant.RotSpeed, 'r');
ylabel('RotSpeed [rpm]');

MyAxes(5) = subplot(n,1,5);
hold on; grid on; box on
plot(FB_Constant.Time,       FB_Constant.GenPwr/1e3, 'r');
ylabel('GenPwr [MW]');

xlabel('time [s]')
linkaxes(MyAxes,'x');    
    

 