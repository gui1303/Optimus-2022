function [Parameter] = NREL5MWDefaultParameter_FBSWE(Parameter)

%% Pitch Controller

% calculated for D_d = [0.700000 ]
% and omega_0_d = [0.500000 ]
% at v_0   = [12.000000 16.000000 20.000000 24.000000 ] 
theta      = [0.095867 0.185201 0.253213 0.321743  ];
kp         = [0.004985 0.004602 0.005397 0.006839 ];
Ti         = [22.471531 32.922070 64.367373 94.532398 ];
ki         = [0.000222 0.000140 0.000084 0.000072 ];
% ---
Parameter.CPC.GS.theta                  = theta;                            % [rad]
Parameter.CPC.GS.kp                     = kp;                               % [s]
Parameter.CPC.GS.Ti                     = Ti;                               % [s] 
Parameter.CPC.GS.ki                     = ki;

Parameter.CPC.Omega_g_rated             = 47.1239;               % [rad/s]
Parameter.CPC.theta_max                 = deg2rad(90);                      % [rad]
Parameter.CPC.theta_min                 = deg2rad(0);                       % [rad]

%% Torque Controller
P_el_rated                              = 17e6;                              % [W]

Parameter.VSC.k                         = 258.52332;                           % [Nm/(rad/s)^2]
Parameter.VSC.Mode                      = 2;                                % 1: ISC, constant power in Region 3; 2: ISC, constant torque in Region 3 
Parameter.VSC.M_g_rated                 = P_el_rated/Parameter.Generator.eta_el/Parameter.CPC.Omega_g_rated;  % [Nm] 
Parameter.VSC.P_a_rated                 = P_el_rated/Parameter.Generator.eta_el;  % [W]


Parameter.VSC.Omega_g_1d5               = 17.80235; % [rad/s];  
Parameter.VSC.kp                        = 141409.989062;
Parameter.VSC.Ti                        = 3.064252 ; 
Parameter.VSC.ki                        = 46148.284973; 

Parameter.VSC.M_g_max                   = Parameter.VSC.M_g_rated*1.1;      % [Nm] 

%% Set-Point-Fading
Parameter.VSC.Delta_Omega_g             = 0.10*Parameter.CPC.Omega_g_rated; % [rad/s]   % over-/under-speed limit for setpoint-fading, first guess
Parameter.VSC.Delta_theta               = deg2rad(20);                      % [rad]     % change of pitch angle at which under-speed limit should be reached, brute-force-optimized
Parameter.VSC.Delta_P                   = 17e6;                          	% [W]       % change of power at which over-speed limit should be reached, first guess
Parameter.Filter.LowPass2.f_cutoff     	= 0.1;                              % [Hz]      % cut-off-frequency for low pass filter for setpoint-fading, first guess


%% Filter Generator Speed
Parameter.Filter.LowPass.Enable         = 1;
Parameter.Filter.LowPass.f_cutoff       = 2;                                % [Hz]

Parameter.Filter.NotchFilter.Enable   	= 1;
Parameter.Filter.NotchFilter.f        	= 1.66;                             % [Hz]
Parameter.Filter.NotchFilter.BW      	= 0.40;                             % [Hz]
Parameter.Filter.NotchFilter.D       	= 0.01;                             % [-]  
   
end