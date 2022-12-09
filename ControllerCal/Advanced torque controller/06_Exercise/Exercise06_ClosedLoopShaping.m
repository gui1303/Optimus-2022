% -----------------------------
% Script: Closed Loop Shaping of Pitch Controller
% Exercise 03 of Master Course 
% "Controller Design for Wind Turbines and Wind Farms"
% ------------
% Task:
% - Get operation point
% - Linearize at each operation point
% - Determine theta, kp and Ti for each operation point
% - Copy the output into NREL5MWDefaultParameter_FBNREL_Ex3.m
% ------------
% History:
% v01:	David Schlipf on 06-Oct-2019
% ----------------------------------

clearvars;close all;clc;

%% Design
OPs         =11;
D_d         = 0.7;
omega_d     = 0.5;

%% Default Parameter Turbine and Controller
Parameter                       = NREL5MWDefaultParameter_SLOW2DOF;
Parameter.VSC.P_a_rated       	= 5e6/Parameter.Generator.eta_el;   % [W]
SteadyStates                    = load('SteadyStatesNREL5MW_FBSWE_SLOW','v_0','Omega','theta'); 
Parameter.Turbine.i = 1/57.35;
Parameter.Turbine.R = 127.5;
Parameter.Turbine.HubHeight = 150;
Parameter.Generator.M_g_dot_max = 150000; % find a better value
Parameter.Generator.eta_el = 0.955;
Parameter.CPC.Omega_g_rated  = 47.1239;
Parameter.VSC.k = 173.24;
Parameter.VSC.M_g_rated = 0.3778e6;
Parameter.VSC.P_a_rated = 17e6;
Parameter.VSC.Omega_g_1d5 = 81.2625/1.5; % improve
Parameter.VSC.M_g_max = 1.1*0.3778e6;
Parameter.VSC.v_rated = 11;
Parameter.VSC.Omega_g_1d5  = 54.175/1.5;  % improve

%% loop over operation points
nOP     = length(OPs);
kp      = NaN(1,nOP);
Ti      = NaN(1,nOP);
theta   = NaN(1,nOP);

for iOP=1:nOP  
    
    % Get operation point
    OP = OPs(iOP);
    
    % Linearize at the operation point
    index = find(abs(SteadyStates.v_0 - OP) < 0.001);
    v_0_OP = OPs(iOP);
    Omega_OP = 0.8221;
    theta_OP = 0;
    [a,b,c,d] = LinearizeSLOW1DOF_TC(Omega_OP,v_0_OP,Parameter);
    b1 = b(1);
    b2 = b(2);
    
    % Determine theta, kp and Ti for the operation point
    kp(iOP) = ((-2*D_d*omega_d+a)/(b1*c));
    ki = -omega_d^2/(b1*c)
    Ti(iOP) = kp(iOP)/ki;
    theta(iOP) = theta_OP;
    
    
    
end

fprintf('Parameter.CPC.GS.theta                  = [%s];\n',sprintf('%f ',theta));
fprintf('Parameter.CPC.GS.kp                     = [%s];\n',sprintf('%f ',kp));
fprintf('Parameter.CPC.GS.Ti                     = [%s];\n',sprintf('%f ',Ti));  
 
