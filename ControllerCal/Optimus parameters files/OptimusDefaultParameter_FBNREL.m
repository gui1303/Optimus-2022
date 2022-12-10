function [Parameter] = NREL5MWDefaultParameter_FBNREL(Parameter)

%% FBSWE Pitch Controller
% calculated for D_d = [0.700000 ]
% and omega_0_d = [0.500000 ]
% at v_0   = [12.000000 16.000000 20.000000 24.000000 ] 
theta      = [0.066817 0.210442 0.304976 0.389997 ];
kp         = [0.066817 0.210442 0.304976 0.389997 ];
Ti         = [3.589214 6.374806 10.558934 15.856681 ];
% ---
Parameter.CPC.GS.theta                  = theta;                            % [rad]
Parameter.CPC.GS.kp                     = kp;                               % [s]
Parameter.CPC.GS.Ti                     = Ti;                               % [s] 

Parameter.CPC.Omega_g_rated             = 47.1239;               % [rad/s]
Parameter.CPC.theta_max                 = deg2rad(90);                      % [rad]
Parameter.CPC.theta_min                 = deg2rad(0);                       % [rad]

%% FBNREL Torque Controller
P_el_rated                              = 17e6;                              % [W]

Parameter.VSC.k                         = 173.24;                           % [Nm/(rad/s)^2]
Parameter.VSC.theta_fine                = deg2rad(1);                       % [rad]      
Parameter.VSC.Mode                      = 1;                                % 1: ISC, constant power in Region 3; 2: ISC, constant torque in Region 3 
Parameter.VSC.M_g_rated                 = P_el_rated/Parameter.Generator.eta_el/Parameter.CPC.Omega_g_rated;  % [Nm] 
Parameter.VSC.P_a_rated                 = P_el_rated/Parameter.Generator.eta_el;  % [W]

% region limits & region parameters based on Jonkman 2009
Parameter.VSC.Omega_g_1To1_5            = 0.59*rpm2radPs(670);                   % [rad/s]
Parameter.VSC.Omega_g_1_5To2            = 0.59*rpm2radPs(871);                   % [rad/s]
Parameter.VSC.Omega_g_2To2_5            = 0.59*rpm2radPs(1150.9);              	% [rad/s]
Parameter.VSC.Omega_g_2_5To3            = 0.59*Parameter.CPC.Omega_g_rated;      % [rad/s]

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

%% Tower Damper
Parameter.TD.gain                       = 0.0436;
Parameter.TD.Power                      = [0 0.8 1 2]*P_el_rated;           % [W]
Parameter.TD.Value                      = [0 0 1 1];

%% Filter Generator Speed
Parameter.Filter.LowPass.Enable         = 1;
Parameter.Filter.LowPass.f_cutoff       = 2;                                % [Hz]

Parameter.Filter.NotchFilter.Enable   	= 1;
Parameter.Filter.NotchFilter.f        	= 1.66;                             % [Hz]
Parameter.Filter.NotchFilter.BW      	= 0.40;                             % [Hz]
Parameter.Filter.NotchFilter.D       	= 0.01;                             % [-]  
                             
Parameter.Filter.BandPass.Enable        = 1;
Parameter.Filter.BandPass.f             = 0.300;                            % [Hz]
Parameter.Filter.BandPass.BW            = 0.10;                             % [Hz]

Parameter.Filter.LowPass2.Enable       	= 1;
Parameter.Filter.LowPass2.f_cutoff     	= 0.1;                              % [Hz]

end