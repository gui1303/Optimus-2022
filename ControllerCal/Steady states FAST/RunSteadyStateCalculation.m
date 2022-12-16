% Steady State Calculation for the IEA 15 MW monopile
% Authors:
% David Schlipf, Mehdy Gooya, Feng Guo
% Copyright (c) 2022 Flensburg University of Applied Sciences, WETI

%% Setup
clearvars;
close all;
clc;
addpath('D:\GitHub\Optimus-2022\MatlabFunctions')
addpath('D:\GitHub\Optimus-2022\OpenFAST modified controller')
addpath('..\..\OpenFAST modified controller');
addpath('..\..\MatlabFunctions');

% Parameters (can be adjusted, but will provide different results)
SaveFlag                = true;    % [true/false]  flag to overwrite SteadyStatesIEA15MW_Monopile_ROSCO_FAST
PlotTimeSignalsFlag     = false;     % [true/false]  flag to plot time results (might be too much for a lot of simulations)    
AdjustSteadyStatesFlag  = true;     % [true/false]  flag to load steady states (set false in the first iteration)
HWindSpeed_vec          = 4:0.1:30;   % [m/s]         range of wind speeds (operation points)
n_Rotation              = 6;        % [-]           number of rotations considered
SteadyStateFile         = 'SteadyStatesOptimus.mat';
Info                    = 'Created by RunSteadyStateCalculation.m of example IEA15MW_09.'; 

% Files (should not be be changed)
FASTexeFile             = 'openfast_x64.exe';
FASTmapFile             = 'MAP_x64.lib';
SimulationName          = 'IEA-15-255-RWT-UMaineSemi';
EDFile                  = 'IEA-15-255-RWT-UMaineSemi_ElastoDyn.dat';
InflowFile              = 'IEA-15-255-RWT_UMaineSemi_InflowFile.dat';
if ~exist('SimulationResults','dir')
    mkdir SimulationResults
end

%% Preprocessing
% change to  WindType 1, simulation length and disable LAC
ManipulateTXTFile(InflowFile,...
                    '4                      WindType',...
                    '1                      WindType');
ManipulateTXTFile([SimulationName,'.fst'],'580   TMax','120   TMax');
%ManipulateTXTFile('ROSCO_v2d6.IN','1 ! FlagLAC','0 ! FlagLAC');

if AdjustSteadyStatesFlag
	load(SteadyStateFile,'v_0','theta','Omega','x_T');
end

%% Processing: run simulations

% Copy the adequate OpenFAST version to the example folder
copyfile(['..\..\OpenFAST modified controller\',FASTexeFile],FASTexeFile)
copyfile(['..\..\OpenFAST modified controller\',FASTmapFile],FASTmapFile)

% Simulate for different steady wind speeds and calculate steady state values
n_HWindSpeed     	= length (HWindSpeed_vec);

for i_HWindSpeed    = 1:n_HWindSpeed
    
    % define FASTresultFile
    HWindSpeed      = HWindSpeed_vec(i_HWindSpeed);
	FASTresultFile 	= ['SimulationResults\',SimulationName,'_',...
     	replace(num2str(HWindSpeed,'%04.1f'),'.','d'),'.outb'];
    
    if ~exist(FASTresultFile,'file')  % run only if simulation does not exist already  
        
        % Adjust the InflowWind file
        ManipulateTXTFile(InflowFile,...
                        '10.0                   HWindSpeed',...
                        [num2str(HWindSpeed,'%4.1f'),' HWindSpeed']);

        % Adjust the ElastoDyn File
        if AdjustSteadyStatesFlag
            MyBlPitch   = num2str(0.0001 + rad2deg  (interp1(v_0,theta,HWindSpeed)),'%5.4f');
            MyRotSpeed  = num2str(0.0001 + radPs2rpm(interp1(v_0,Omega,HWindSpeed)),'%5.4f');
            MyTTDspFA   = num2str(0.0001  +          (interp1(v_0,x_T  ,HWindSpeed)),'%5.4f');
        else
            MyBlPitch   = num2str(0.9875);
            MyRotSpeed  = num2str(7.85);
            MyTTDspFA   = num2str(-0.2754);
        end               
        ManipulateTXTFile(EDFile,'MyBlPitch', MyBlPitch);
        ManipulateTXTFile(EDFile,'MyRotSpeed',MyRotSpeed);
        ManipulateTXTFile(EDFile,'MyTTDspFA', MyTTDspFA);                 

        % Run FB 
        dos([FASTexeFile,' ',SimulationName,'.fst']);
        movefile([SimulationName,'.outb'],FASTresultFile)

        % Reset the InflowWind file again
        ManipulateTXTFile(InflowFile,...                    
                        [num2str(HWindSpeed,'%4.1f'),' HWindSpeed'],...
                        '10.0                   HWindSpeed');
                    
        % Reset the ElastoDyn file again
        ManipulateTXTFile(EDFile,MyBlPitch,'MyBlPitch');
        ManipulateTXTFile(EDFile,MyRotSpeed,'MyRotSpeed');
        ManipulateTXTFile(EDFile,MyTTDspFA,'MyTTDspFA');
    end
end

% Reset WindType, simulation length and LAC
ManipulateTXTFile('IEA-15-255-RWT_UMaineSemi_InflowFile.dat',...
                    '1                      WindType',...
                    '4                      WindType');
ManipulateTXTFile([SimulationName,'.fst'],'120   TMax','580   TMax');
%ManipulateTXTFile('ROSCO_v2d6.IN','0 ! FlagLAC','1 ! FlagLAC');

% Clean up
delete(FASTexeFile)
delete(FASTmapFile)

%% Postprocessing: generate steady state file

% allocation SLOW style steady states
v_0             = NaN(1,n_HWindSpeed);
Omega           = NaN(1,n_HWindSpeed);
theta           = NaN(1,n_HWindSpeed);
x_T             = NaN(1,n_HWindSpeed);
M_g             = NaN(1,n_HWindSpeed);

% allocation of evaluation values
STD_Omega       = NaN(1,n_HWindSpeed);
STD_x_T         = NaN(1,n_HWindSpeed);
RotSpeedCell    = cell(1,n_HWindSpeed);
TTDspFACell     = cell(1,n_HWindSpeed);

% loop over all wind speeds
for i_HWindSpeed = 1:n_HWindSpeed
    
    % Load data
    HWindSpeed              = HWindSpeed_vec(i_HWindSpeed);
    FASTresultFile          = ['SimulationResults\',SimulationName,'_',...
        replace(num2str(HWindSpeed,'%04.1f'),'.','d'),'.outb'];    
    FB                      = ReadFASTbinaryIntoStruct(FASTresultFile);

    % estimation of rotor speed, one iteration
    Omega_temp              = rpm2radPs(FB.RotSpeed(end));
    T_Rotation              = 2*pi/Omega_temp;
    Considered              = FB.Time>FB.Time(end)-n_Rotation*T_Rotation;
    Omega_temp              = rpm2radPs(mean(FB.RotSpeed(Considered)));
    T_Rotation              = 2*pi/Omega_temp;
    Considered              = FB.Time>FB.Time(end)-n_Rotation*T_Rotation;
    
    % calculate steady state values as mean over last n_Rotation       
    v_0(i_HWindSpeed)       = HWindSpeed;
    theta(i_HWindSpeed)     = deg2rad  (mean(FB.BldPitch1(Considered)));
    M_g(i_HWindSpeed)       = 1e3*     (mean(FB.GenTq    (Considered)));
    Omega(i_HWindSpeed)     = rpm2radPs(mean(FB.RotSpeed (Considered)));
    x_T(i_HWindSpeed)       =          (mean(FB.TTDspFA  (Considered)));    
    
	% Calculate standard deviation
    STD_Omega(i_HWindSpeed)	= rpm2radPs(std(FB.RotSpeed  (Considered)));
    STD_x_T(i_HWindSpeed)	=          (std(FB.TTDspFA   (Considered)));
    
    % Store time signals, if requested
    if PlotTimeSignalsFlag
        RotSpeedCell{i_HWindSpeed}  = FB.RotSpeed;
        TTDspFACell {i_HWindSpeed}  = FB.TTDspFA;
    end

end

% Save if requested
if SaveFlag
    save(SteadyStateFile,'v_0','theta','Omega','x_T','M_g','Info');
end

% Plot config
MyMarkerSize    = 20;
ScreenSize      = get(0,'ScreenSize');

% Plot std
figure('Name','Steady States','position',[.1 .1 .8 .8].*ScreenSize([3,4,3,4]))
n       = 2;

MyAxes(1) = subplot(n,1,1);
hold on; grid on; box on
plot(v_0,radPs2rpm(STD_Omega),'.-','MarkerSize',MyMarkerSize);
ylabel({'STD(RotSpeed)';'[rpm]'});

MyAxes(2) = subplot(n,1,2);
hold on; grid on; box on
plot(v_0,STD_x_T,'.-','MarkerSize',MyMarkerSize);
ylabel({'STD(TTDspFA)';'[m]'});

xlabel('wind speed [m/s]')
linkaxes(MyAxes,'x');

% Plot steady states
clear MyAxes
figure('Name','Steady States','position',[.1 .1 .8 .8].*ScreenSize([3,4,3,4]))
n               = 4;

MyAxes(1) = subplot(n,1,1);
hold on; grid on; box on
plot(v_0,rad2deg(theta),'.-','MarkerSize',MyMarkerSize);
ylabel({'BldPitch1'; '[deg]'});

MyAxes(2) = subplot(n,1,2);
hold on; grid on; box on
plot(v_0,M_g,'.-','MarkerSize',MyMarkerSize);
ylabel({'GenTq';'[MNm]'});

MyAxes(3) = subplot(n,1,3);
hold on; grid on; box on
plot(v_0,radPs2rpm(Omega),'.-','MarkerSize',MyMarkerSize);
ylabel({'RotSpeed';'[rpm]'});

MyAxes(4) = subplot(n,1,4);
hold on; grid on; box on
plot(v_0,x_T,'.-','MarkerSize',MyMarkerSize);
ylabel({'TTDspFA';'[m]'});

xlabel('wind speed [m/s]')
linkaxes(MyAxes,'x');

% Plot time signals, if requested

if PlotTimeSignalsFlag
    clear MyAxes
    figure('Name','Steady States','position',[.1 .1 .8 .8].*ScreenSize([3,4,3,4]))
    n       = 2;
    
    MyAxes(1) = subplot(n,1,1);
    hold on; grid on; box on
    plot(FB.Time,cell2mat(RotSpeedCell));
    ylabel({'RotSpeed';'[rpm]'});

    MyAxes(2) = subplot(n,1,2);
    hold on; grid on; box on
    plot(FB.Time,cell2mat(TTDspFACell));
    ylabel({'TTDspFA';'[m]'});

    xlabel('time [s]')
    linkaxes(MyAxes,'x');
end
