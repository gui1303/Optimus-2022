% -----------------------------
% Function: provides parameter for NREL5MW SLOW wind turbine model
% ------------
% Input:
% none
% ------------
% Output:
% - Parameter   struct of Parameters
% ------------
% History: 
% v01:	David Schlipf on 29-Sep-2019
% ----------------------------------
function Parameter = NREL5MWDefaultParameter_SLOW1DOF

%% General          
Parameter.General.rho               = 1.225;         	% [kg/m^3]  air density

%% Turbine
Parameter.Turbine.i             	= 1/97;             % [-]       gearbox ratio
Parameter.Turbine.R              	= 126/2;            % [m]       Rotor radius
Parameter.Turbine.SS             	= load('PowerAndThrustCoefficientsNREL5MW','c_P','theta','lambda'); % load Power coefficient look-up-table
J_G                               	= 534.116;          % [kgm^2]	generator inertia about high-speed shaft
J_R                                	= 3.8759e+007;      % [kgm^2]	rotor inertia about low-speed shaft
Parameter.Turbine.J                	= J_R+J_G/Parameter.Turbine.i^2; % [kgm^2] sum of moments of inertia about low-speed shaft

%% Generator
Parameter.Generator.eta_el      	= 0.944;            % [-]       Generator efficency
Parameter.Generator.M_g_dot_max     = 15e3;             % [Nm/s]    Maximum Torque rate

end