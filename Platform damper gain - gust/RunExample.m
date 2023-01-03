% Authors:
% David Schlipf, Feng Guo

% Modified to load steady states and optimize platform gain value to
% minimize standard deviation
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


% Number of gains to be tried out
vGains = -6:0.25:-4;
nGains = length(vGains);
vSTD   = NaN(1,nGains);
vMean  = NaN(1,nGains);

if ~exist('PtfmPitchInfo.mat')

    for iGain = 1:nGains

        % Adjust ROSCO
        Gain = num2str(vGains(iGain));
        ManipulateTXTFile('ROSCO_v2d6.IN','MyGain         ! Fl_Kp' ,[Gain,'       ! Fl_Kp']); 

        % Load steady states and ajust initial conditions
        HWindSpeed              = 19; % Has to match gust wnd file (+1)
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
        FASTresultFile      = 'IEA-15-255-RWT-UMaineSemi.outb';    
        dos([FASTexeFile,' ',SimulationName,'.fst']);
        FB_Gust             = ReadFASTbinaryIntoStruct(FASTresultFile);


        delete(FASTexeFile)
        delete(FASTmapFile)
        delete(FASTresultFile)

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

        % Reset ROSCO
        ManipulateTXTFile('ROSCO_v2d6.IN',[Gain,'       ! Fl_Kp'],'MyGain         ! Fl_Kp'); 

        % Calculate standard deviation for each gain
        vSTD(iGain)  = std(FB_Gust.PtfmPitch);
        vMean(iGain) = mean(FB_Gust.PtfmPitch);

    end
    info = num2str(vGains);
    save('PtfmPitchInfo.mat','vSTD','vMean','info');
    
else
    load('PtfmPitchInfo.mat','vSTD','vMean');
end


 