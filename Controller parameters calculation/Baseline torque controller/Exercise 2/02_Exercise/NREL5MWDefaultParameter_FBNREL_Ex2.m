% -----------------------------
% Function: should add parameters for NREL5MW Baseline Torque Controller
% Exercise 02 of Master Course 
% "Controller Design for Wind Turbines and Wind Farms"
% ------------
% Input:
% - Parameter   struct of Parameters
% ------------
% Output:
% - Parameter   struct of Parameters
% ------------
% History: 
% v02:	David Schlipf on 26-Jul-2021: small update
% v01:	David Schlipf on 29-Sep-2019
% ----------------------------------
function Parameter = NREL5MWDefaultParameter_FBNREL_Ex2(Parameter)

%% FBNREL Torque Controller
Omega_g_rated                           = rpm2radPs(12.1*97);               % [rad/s]
P_el_rated                              = 5e6;                              % [W]
lambda_opt                              = 7.55;                              % [-]
theta_opt                               = deg2rad(0);                              % [deg]
% Use interp2 for CPopt
c_P_opt                                 = interp2(Parameter.Turbine.SS.theta,Parameter.Turbine.SS.lambda,Parameter.Turbine.SS.c_P,theta_opt,lambda_opt);
rho                                     = Parameter.General.rho;
R                                       = Parameter.Turbine.R;
i                                       = Parameter.Turbine.i;
Parameter.VSC.k                         = (1/2)*rho*pi*R^5*(c_P_opt/lambda_opt^3)*i^3;                              % [Nm/(rad/s)^2]
Parameter.VSC.theta_fine                = deg2rad(1);                       % [rad]      
Parameter.VSC.Mode                      = 1;                                % 1: ISC, constant power in Region 3; 2: ISC, constant torque in Region 3 
Parameter.VSC.P_a_rated                 = P_el_rated/Parameter.Generator.eta_el;  % [W] aerodynamic power
Parameter.VSC.M_g_rated                 = P_el_rated/(Omega_g_rated*Parameter.Generator.eta_el);                              % [Nm] 

% region limits & region parameters based on Jonkman 2009
Parameter.VSC.Omega_g_1To1_5            = rpm2radPs(670);                   % [rad/s]
Parameter.VSC.Omega_g_1_5To2            = rpm2radPs(871);                   % [rad/s]
Parameter.VSC.Omega_g_2To2_5            = rpm2radPs(1150.9);              	% [rad/s]
Parameter.VSC.Omega_g_2_5To3            = Omega_g_rated;                    % [rad/s]

% Region 1_5: M_g = a * Omega_g + b: 
% 1.Eq: 0                   = a * Omega_g_1To1_5 + b 
% 2.Eq: k*Omega_g_1_5To2^2  = a * Omega_g_1_5To2 + b
Parameter.VSC.a_1_5                     = 918.1;
Parameter.VSC.b_1_5                     = -64417;

% Region 2_5: M_g = a * Omega_g + b: 
% 1.Eq: M_g_rated           = a * Omega_g_2_5To3   	+ b 
% 2.Eq: k*Omega_g_2To2_5^2  = a * Omega_g_2To2_5    + b 
Parameter.VSC.a_2_5                     = 3977.9;
Parameter.VSC.b_2_5                     = -445606;

%% Calculating best TSR
[CP_max,index] = max(Parameter.Turbine.SS.c_P(:,1));
TSR_max = Parameter.Turbine.SS.lambda(index); 
PowerIncrease = ((CP_max / c_P_opt) - 1) * 100;

end