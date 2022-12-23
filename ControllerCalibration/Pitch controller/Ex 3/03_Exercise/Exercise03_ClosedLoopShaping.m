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
OPs         = [12 14 18 24];
D_d         = 1.2;
omega_d     = 0.15;

%% Default Parameter Turbine and Controller
Parameter                       = OptimusDefaultParameter_SLOW1DOF;
Parameter.VSC.P_a_rated       	= 17e6/Parameter.Generator.eta_el;   % [W]
SteadyStates                    = load('SteadyStatesOptimusFAST','v_0','Omega','theta');                       

%% loop over operation points
nOP     = length(OPs);
kp      = NaN(1,nOP);
Ti      = NaN(1,nOP);
theta   = NaN(1,nOP);

for iOP=1:nOP  
    
    % Get operation point
    OP = OPs(iOP);
    
    % Linearize at each operation point
    index = find(abs(SteadyStates.v_0 - OP) < 0.043);
    v_0_OP = OPs(iOP);
    Omega_OP = SteadyStates.Omega(index);
    theta_OP = SteadyStates.theta(index);
    [a,b,c,d] = LinearizeSLOW1DOF_PC(theta_OP,Omega_OP,v_0_OP,Parameter);
    b1 = b(1);
    b2 = b(2);
    
    % Determine theta, kp and Ti for each operation point
    kp(iOP) = ((-2*D_d*omega_d+a)/(b1*c));
    ki(iOP) = -omega_d^2/(b1*c);
    Ti(iOP) = kp(iOP)/ki(iOP);
    theta(iOP) = theta_OP;
    
    
    
end

fprintf('Parameter.CPC.GS.theta                  = [%s];\n',sprintf('%f ',theta));
fprintf('Parameter.CPC.GS.kp                     = [%s];\n',sprintf('%f ',kp));
fprintf('Parameter.CPC.GS.Ti                     = [%s];\n',sprintf('%f ',Ti));  
fprintf('Parameter.CPC.GS.ki                     = [%s];\n',sprintf('%f ',ki)); 