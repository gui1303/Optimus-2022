% -----------------------------
% Function: should add parameters for NREL5MW Pitch Controller
% Exercise 03 of Master Course 
% "Controller Design for Wind Turbines and Wind Farms"
% ------------
% Input:
% - Parameter   struct of Parameters
% ------------
% Output:
% - Parameter   struct of Parameters
% ------------
% History: 
% v01:	David Schlipf on 06-Oct-2019
% ----------------------------------
function [Parameter] = NREL5MWDefaultParameter_FBNREL_Ex3(Parameter)

%% FBNREL Torque Controller
Omega_g_rated                           = rpm2radPs(12.1*97);               % [rad/s]
P_el_rated                              = 5e6;                              % [W]
lambda_opt                              = 7.55;                             % [-]
c_P_opt                                 = interp2(Parameter.Turbine.SS.theta,Parameter.Turbine.SS.lambda,Parameter.Turbine.SS.c_P,0,lambda_opt);
rho                                     = Parameter.General.rho;
R                                       = Parameter.Turbine.R;
i                                       = Parameter.Turbine.i;  
Parameter.VSC.k                         = 1/2*rho*pi*R^5*c_P_opt/lambda_opt^3*i^3;  % [Nm/(rad/s)^2]
Parameter.VSC.theta_fine                = deg2rad(1);                       % [rad]      
Parameter.VSC.Mode                      = 1;                                % 1: ISC, constant power in Region 3; 2: ISC, constant torque in Region 3 
Parameter.VSC.M_g_rated                 = P_el_rated/Parameter.Generator.eta_el/Omega_g_rated;  % [Nm] 
Parameter.VSC.P_a_rated                 = P_el_rated/Parameter.Generator.eta_el;  % [W]

% region limits & region parameters based on Jonkman 2009
Parameter.VSC.Omega_g_1To1_5            = rpm2radPs(670);                   % [rad/s]
Parameter.VSC.Omega_g_1_5To2            = rpm2radPs(871);                   % [rad/s]
Parameter.VSC.Omega_g_2To2_5            = rpm2radPs(1150.9);              	% [rad/s]
Parameter.VSC.Omega_g_2_5To3            = Omega_g_rated;                    % [rad/s]

% Region 1_5: M_g = a * Omega_g + b: 
% 1.Eq: 0                   = a * Omega_g_1To1_5 + b 
% 2.Eq: k*Omega_g_1_5To2^2  = a * Omega_g_1_5To2 + b
Parameter.VSC.a_1_5                     = Parameter.VSC.k*Parameter.VSC.Omega_g_1_5To2^2/(Parameter.VSC.Omega_g_1_5To2-Parameter.VSC.Omega_g_1To1_5);
Parameter.VSC.b_1_5                     = -Parameter.VSC.a_1_5*Parameter.VSC.Omega_g_1To1_5;

% Region 2_5: M_g = a * Omega_g + b: 
% 1.Eq: M_g_rated           = a * Omega_g_2_5To3   	+ b 
% 2.Eq: k*Omega_g_2To2_5^2  = a * Omega_g_2To2_5    + b
Parameter.VSC.a_2_5                     = (Parameter.VSC.M_g_rated-Parameter.VSC.k*Parameter.VSC.Omega_g_2To2_5^2)/(Parameter.VSC.Omega_g_2_5To3-Parameter.VSC.Omega_g_2To2_5);
Parameter.VSC.b_2_5                     = Parameter.VSC.M_g_rated-Parameter.VSC.a_2_5*Parameter.VSC.Omega_g_2_5To3;


%% FBSWE Pitch Controller
Parameter.CPC.GS.theta                  = [0.066817 0.210442 0.304976 0.389997];
Parameter.CPC.GS.kp                     = [0.020485 0.008687 0.006984 0.007214];
Parameter.CPC.GS.Ti                     = [2.702810 3.167971 3.791372 4.132752];

Parameter.CPC.Omega_g_rated             = Omega_g_rated;                    % [rad/s]
Parameter.CPC.theta_max                 = deg2rad(90);                      % [rad]
Parameter.CPC.theta_min                 = deg2rad(0);                       % [rad]


end