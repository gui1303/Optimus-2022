function [Parameter] = NREL5MWDefaultParameter_FBSWE_Ex6_TPI(Parameter)

%% FBSWE Pitch Controller

Parameter.CPC.Omega_g_rated             = rpm2radPs(12.1*97);               % [rad/s]

%% FBSWE Torque Controller
P_el_rated                              = 5e6;                              % [W]

Parameter.VSC.k                         = 2.3689;                           % [Nm/(rad/s)^2]
Parameter.VSC.Mode                      = 1;                                % 1: ISC, constant power in Region 3; 2: ISC, constant torque in Region 3 
Parameter.VSC.M_g_rated                 = P_el_rated/Parameter.Generator.eta_el/Parameter.CPC.Omega_g_rated;  % [Nm] 
Parameter.VSC.P_a_rated                 = P_el_rated/Parameter.Generator.eta_el;  % [W]
Parameter.VSC.Omega_g_1d5               = rpm2radPs(8)/Parameter.Turbine.i; % [rad/s];  
Parameter.VSC.M_g_max                   = Parameter.VSC.M_g_rated*1.1;      % [Nm] 

Parameter.VSC.kp                        = 3540.733821;                                % [Nm/(rad/s)]  % Please adjust !!!
Parameter.VSC.Ti                        = 3.043521;                                % [s]           % Please adjust !!!

%% Filter Generator Speed
Parameter.Filter.LowPass.Enable         = 1;
Parameter.Filter.LowPass.f_cutoff       = 2;                                % [Hz]

Parameter.Filter.NotchFilter.Enable   	= 1;
Parameter.Filter.NotchFilter.f        	= 1.66;                             % [Hz]
Parameter.Filter.NotchFilter.BW      	= 0.40;                             % [Hz]
Parameter.Filter.NotchFilter.D       	= 0.01;                             % [-]  
   
end