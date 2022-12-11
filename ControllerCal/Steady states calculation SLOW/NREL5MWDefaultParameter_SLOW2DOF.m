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
% v01:	David Schlipf on 20-Oct-2019
% ----------------------------------
function [Parameter] = NREL5MWDefaultParameter_SLOW2DOF

%% General          
Parameter.General.rho               = 1.225;         	% [kg/m^3]  air density
%--------------------------------------------------------------------------
%% Turbine
Parameter.Turbine.SS             	= load('PowerAndThrustCoefficientsNREL5MW','c_P','c_T','theta','lambda');

Parameter.Turbine.i               	= 1/57.35;             % [-]       gearbox ratio
Parameter.Turbine.R              	= 255/2;            % [m]       Rotor radius

% Generator
Parameter.Generator.eta_el      	= 0.955;            % [-]       Generator efficiency
Parameter.Generator.M_g_dot_max     = 0.3778e6/4;             % [Nm/s]    maximum torque rate, 1/4 of rated

% drive-train dynamics
J_G                               	= 3*534.116;          % [kgm^2]	Generator Inertia About High-Speed Shaft
J_R                                	= 3*3.8759e+007;      % [kgm^2]	Rotor Inertia About High-Speed Shaft
Parameter.Turbine.J                	= J_R+J_G/Parameter.Turbine.i^2;

% fore-aft tower dynamics  
d_s                                 = 0.01;             % [-]       Structural Damping ratio from NRELOffshrBsline5MW_Tower_Onshore.dat
f_0TwFADOF1                         = 0.324;            % [Hz]      first tower fore-aft eigenfrequency
Parameter.Turbine.x_T0           	= 2*(-0.0140);          % [m]       tower top deflection without wind
Parameter.Turbine.k_Te            	= 1.81e+06;         % [kg/s^2]  tower equivalent bending stiffness 
Parameter.Turbine.m_Te          	= Parameter.Turbine.k_Te/(f_0TwFADOF1*2*pi)^2;  	% [kg]      tower equivalent modal mass
Parameter.Turbine.c_Te            	= 2*Parameter.Turbine.m_Te*d_s*f_0TwFADOF1*2*pi;    % [kg/s]	tower equivalent structual damping
Parameter.Turbine.HubHeight       	= 150;               % [m]       hub height

end