% Copy of RunExample from the IEW15MW_02 from the 
% https://github.com/fengguoFUAS/Baseline-Lidar-assisted-Controller/commit/0d83d87b6c7b03013fc2ad6e12c4fe993419dbd2)
% Modified to a simpler version without lidar and FF control.

%% Setup
clearvars; 
close all; 
clc;
addpath('..\MatlabFunctions')

%% Run FB without platform damper

% Disable the damper. Should have some kind of test here in case the damper
% is already disabled it will not change.
ManipulateTXTFile('ROSCO_v2d6.IN','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');

% Copy the adequate OpenFAST version to the example folder
FASTexeFile     = 'openfast_x64.exe';
FASTmapFile     = 'MAP_x64.lib';
copyfile(['..\OpenFAST modified controller\',FASTexeFile],FASTexeFile)
copyfile(['..\OpenFAST modified controller\',FASTmapFile],FASTmapFile)

% Run the simulation
SimulationName  = 'IEA-15-240-RWT-UMaineSemiNoDamper';
dos([FASTexeFile,' ',SimulationName,'.fst']);                       % run OpenFAST
movefile([SimulationName,'.outb'],[SimulationName,'_FB.outb'])      % store results

% read in data
FB_no_damper     = ReadFASTbinaryIntoStruct([SimulationName,'_FB.outb']);

% Clean up
delete(FASTexeFile)
delete(FASTmapFile)

%% Run FB with damper

% Enable the damper
ManipulateTXTFile('ROSCO_v2d6.IN','0                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}','3                   ! Fl_Mode           - Floating specific feedback mode {0: no nacelle velocity feedback, 1: feed back translational velocity, 2: feed back rotational veloicty}');

% Copy again the adequate OpenFAST version to the example folder
FASTexeFile     = 'openfast_x64.exe';
FASTmapFile     = 'MAP_x64.lib';
copyfile(['..\OpenFAST modified controller\',FASTexeFile],FASTexeFile)
copyfile(['..\OpenFAST modified controller\',FASTmapFile],FASTmapFile)
% Run the simulation
SimulationName  = 'IEA-15-240-RWT-UMaineSemiWithDamper';
dos([FASTexeFile,' ',SimulationName,'.fst']);                       % run OpenFAST
movefile([SimulationName,'.outb'],[SimulationName,'_FB.outb'])      % store results
% read in data
FB_damper     = ReadFASTbinaryIntoStruct([SimulationName,'_FB.outb']);
% Clean up
delete(FASTexeFile)
delete(FASTmapFile)

%% Comparison


% Plot         
ScreenSize = get(0,'ScreenSize');
figure('Name','Simulation results','position',[.1 .1 .8 .8].*ScreenSize([3,4,3,4]))

MyAxes(1) = subplot(4,1,1);
hold on; grid on; box on
plot(FB_no_damper.Time,       FB_no_damper.Wind1VelX);
plot(FB_damper.Time,       FB_damper.Wind1VelX);
ylabel('[m/s]');

MyAxes(2) = subplot(4,1,2);
hold on; grid on; box on
plot(FB_no_damper.Time,       FB_no_damper.BldPitch1, 'r');
plot(FB_damper.Time,       FB_damper.BldPitch1, 'b' );
ylabel('BldPitch1 [deg]');

MyAxes(3) = subplot(4,1,3);
hold on; grid on; box on
plot(FB_no_damper.Time,       FB_no_damper.RotSpeed, 'r');
plot(FB_damper.Time,       FB_damper.RotSpeed, 'b');
ylabel('RotSpeed [rpm]');
legend ('Platform Damper off', 'Platform damper on')

MyAxes(4) = subplot(4,1,4);
hold on; grid on; box on
plot(FB_no_damper.Time,       FB_no_damper.PtfmPitch, 'r');
plot(FB_damper.Time,       FB_damper.PtfmPitch, 'b');
ylabel('PtfmPitch [deg]');

xlabel('time [s]')
linkaxes(MyAxes,'x');