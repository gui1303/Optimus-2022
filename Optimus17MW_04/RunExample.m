% Authors:
% David Schlipf, Feng Guo

% Modified to load steady states 
% Copyright (c) 2022 Flensburg University of Applied Sciences, WETI

%% Setup
clearvars;
close all;
clc;
addpath('..\MatlabFunctions')
addpath('..\MatlabFunctions\AnalyticlModel')


% Files (should not be be changed)
FASTexeFile         = 'openfast_x64.exe';
FASTmapFile         = 'MAP_x64.lib';
SimulationName      = 'IEA-15-255-RWT-UMaineSemi';

if ~exist('SimulationResultsStep','dir')
    mkdir SimulationResultsStep
end

% Load steady states and ajust initial conditions
HWindSpeed              = 18; % Has to match step wnd file
SteadyStateFile         = 'SteadyStatesOptimusFAST.mat';
EDFile                  = 'IEA-15-255-RWT-UMaineSemi_ElastoDyn.dat';
load(SteadyStateFile,'v_0','theta','Omega','x_T','TdpsSS','PtfSurge','PtfSway','PtfHeave','PtfRoll','PtfPitch','PtfYaw');
MyBlPitch   = num2str(0.0001  + rad2deg  (interp1(v_0,theta,HWindSpeed)),'%5.4f');
MyRotSpeed  = num2str(0.0001  + radPs2rpm(interp1(v_0,Omega,HWindSpeed)),'%5.4f');
MyTTDspFA   = num2str(0.0001  +          (interp1(v_0,x_T  ,HWindSpeed)),'%5.4f');
MyTTDspSS   = num2str(0.0001  +          (interp1(v_0,TdpsSS  ,HWindSpeed)),'%5.4f');
MyPtfmSurge = num2str(0.0001  +          (interp1(v_0,PtfSurge  ,HWindSpeed)),'%5.4f');
MyPtfmSway  = num2str(0.0001  +          (interp1(v_0,PtfSway  ,HWindSpeed)),'%5.4f');
MyPtfmHeave = num2str(0.0001  +          (interp1(v_0,PtfHeave  ,HWindSpeed)),'%5.4f');
MyPtfmRoll  = num2str(0.0001  +          (interp1(v_0,PtfRoll  ,HWindSpeed)),'%5.4f');
MyPtfmPitch = num2str(0.0001  +          (interp1(v_0,PtfPitch  ,HWindSpeed)),'%5.4f');
MyPtfmYaw   = num2str(0.0001  +          (interp1(v_0,PtfYaw  ,HWindSpeed)),'%5.4f');

% Change ICs
ManipulateTXTFile(EDFile,'MyBlPitch', MyBlPitch);
ManipulateTXTFile(EDFile,'MyRotSpeed',MyRotSpeed);
ManipulateTXTFile(EDFile,'MyTTDspFA', MyTTDspFA);  
ManipulateTXTFile(EDFile,'MyTTDspSS', MyTTDspSS); 
ManipulateTXTFile(EDFile,'MyPtfmSurge', MyPtfmSurge); 
ManipulateTXTFile(EDFile,'MyPtfmSway', MyPtfmSway); 
ManipulateTXTFile(EDFile,'MyPtfmHeave', MyPtfmHeave); 
ManipulateTXTFile(EDFile,'MyPtfmRoll', MyPtfmRoll); 
ManipulateTXTFile(EDFile,'MyPtfmPitch', MyPtfmPitch); 
ManipulateTXTFile(EDFile,'MyPtfmYaw', MyPtfmYaw);



%% Processing: run simulations with gust

% Copy the adequate OpenFAST version to the example folder
copyfile(['..\OpenFAST modified controller\',FASTexeFile],FASTexeFile)
copyfile(['..\OpenFAST modified controller\',FASTmapFile],FASTmapFile)


% Run the simulation
FASTresultFile      = 'SimulationResultsStep\URef_18_Step.outb';
if ~exist(FASTresultFile,'file')    
    dos([FASTexeFile,' ',SimulationName,'.fst']);
    movefile([SimulationName,'.outb'],FASTresultFile)
end     

delete(FASTexeFile)
delete(FASTmapFile)

% Reset the ElastoDyn file again
ManipulateTXTFile(EDFile,MyBlPitch,'MyBlPitch');
ManipulateTXTFile(EDFile,MyRotSpeed,'MyRotSpeed');
ManipulateTXTFile(EDFile,MyTTDspFA,'MyTTDspFA');
ManipulateTXTFile(EDFile,MyTTDspFA,'MyTTDspSS');
ManipulateTXTFile(EDFile,MyPtfmSurge,'MyPtfmSurge');
ManipulateTXTFile(EDFile,MyPtfmSway,'MyPtfmSway');
ManipulateTXTFile(EDFile,MyPtfmHeave,'MyPtfmHeave');
ManipulateTXTFile(EDFile,MyPtfmRoll,'MyPtfmRoll');
ManipulateTXTFile(EDFile,MyPtfmPitch,'MyPtfmPitch');
ManipulateTXTFile(EDFile,MyPtfmYaw,'MyPtfmYaw');

% read in data
FB_Gust     = ReadFASTbinaryIntoStruct(FASTresultFile);


% Plot time results
figure('Name','Time results for wind')
hold on; grid on; box on
plot(FB_Gust.Time,FB_Gust.Wind1VelX,'Color',[0.8500 0.3250 0.0980]);
ylabel('Wind velocity [m/s]');
xlabel('time [s]')
xlim([0 630]);
hold off

figure('Name','Time results for rotor speed')
hold on; grid on; box on
plot(FB_Gust.Time,FB_Gust.RotSpeed,'Color',[0.8500 0.3250 0.0980]);
ylabel('Rotor speed');
xlabel('time [s]')
xlim([0 630]);
hold off
    
figure('Name','Time results for platform pitch')
hold on; grid on; box on
plot(FB_Gust.Time,FB_Gust.PtfmPitch,'Color',[0.8500 0.3250 0.0980]);
ylabel('Platform pitch [deg]');
xlabel('time [s]')
xlim([0 630]);
hold off
    
figure('Name','Time results for blade pitch')
hold on; grid on; box on
plot(FB_Gust.Time,FB_Gust.BldPitch1,'Color',[0.8500 0.3250 0.0980]);
ylabel('Blade pitch [deg]');
xlabel('time [s]')
xlim([0 630]);
hold off
    
figure('Name','Time results for electric power')
hold on; grid on; box on
plot(FB_Gust.Time,FB_Gust.GenPwr/1000,'Color',[0.8500 0.3250 0.0980]);
ylabel('Electric power [MW]');
xlabel('time [s]')
xlim([0 630]);
hold off
    
    

 