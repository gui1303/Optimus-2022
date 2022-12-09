% -----------------------------
% Script: Calculates SteadyStates
% Exercise 08 of Master Course 
% "Controller Design for Wind Turbines and Wind Farms"
% ------------
% Task:
% - adjust script and config
% ------------
% History:
% v02:	David Schlipf on 31-Dec-2020
% v01:	David Schlipf on 24-Nov-2019
% ----------------------------------
%% 1. Initialitzation 
clearvars; close all;clc;

%% 2. Config
[v_0,FlagPITorqueControl,Parameter]    = StaticCalculationsConfig_Solutions;

%% 3. Allocation
Omega               = zeros(1,length(v_0));
theta               = zeros(1,length(v_0));
M_g                 = zeros(1,length(v_0));
Omega_dot_Sq        = zeros(1,length(v_0));
exitflag            = zeros(1,length(v_0));


%% Adjust parameters for Optimus 
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
kp =21994.716514;
ki = 3.328078284961565e+03;
%% 4. Loop over wind speeds to determine omega, theta, M_g
for iv_0=1:length(v_0)
    v_0i        	= v_0(iv_0);
    
    %% 4.1 Determine Region
    if FlagPITorqueControl
        if      v_0i < Parameter.VSC.v_1to1d5
            Region = '1';
        elseif  v_0i < Parameter.VSC.v_1d5to2
        	Region = '1.5';
        elseif  v_0i < Parameter.VSC.v_2to2d5    
            Region = '2';
        elseif  v_0i < Parameter.VSC.v_rated 
            Region = '2.5';
        else
            Region = '3';
        end        
      
    else % no PI torque control
        if      v_0i < Parameter.VSC.v_rated
            Region = 'StateFeedback';
        else
            Region = '3';
        end
    end

    %% 4.2 Determine Static Values
    switch Region %
        
        case '1' %  % Determin Omega in Region 1, where theta and M_g are fixed
            
            theta_min   = Parameter.CPC.theta_min;
            Omega_min   = 0;
            Omega_max   = Parameter.VSC.Omega_g_1d5*Parameter.Turbine.i;
            M_g_min     = 0;

            [Omega(iv_0),Omega_dot_Sq(iv_0),exitflag(iv_0)] = ...         
                fminbnd(@(Omega) (OmegaDot(Omega,theta_min,M_g_min,v_0i,Parameter))^2,...
                Omega_min,Omega_max,optimset('Display','none'));

            theta(iv_0) = theta_min;
            M_g(iv_0)   = M_g_min;
            
        case '1.5' % Determin M_g in Region 1.5, where Omega and theta are fixed

            theta_min   = Parameter.CPC.theta_min;
            Omega_1d5   = Parameter.VSC.Omega_g_1d5*Parameter.Turbine.i;
            M_g_min     = 0;
            M_g_max     = Parameter.VSC.M_g_max;

            [M_g(iv_0),Omega_dot_Sq(iv_0),exitflag(iv_0)] = ...         
                fminbnd(@(M_g) (OmegaDot(Omega_1d5,theta_min,M_g,v_0i,Parameter))^2,...
                M_g_min,M_g_max,optimset('Display','none'));

            theta(iv_0) = theta_min;            
            Omega(iv_0) = Omega_1d5;  
            
        case {'2','StateFeedback'} % Determin Omega and M_g in Region 2 (or 1-2.5 for state feedback), where theta is fixed 
            
            Omega_min   = rpm2radPs(7); % to avoid stall solution
            Omega_max   = Parameter.CPC.Omega_g_rated*Parameter.Turbine.i;
            theta_i     = Parameter.CPC.theta_min;

            [Omega_i,Omega_dot_Sq(iv_0),exitflag(iv_0)] = ...
                fminbnd(@(Omega) (OmegaDot(Omega,theta_i,v_0i,Parameter))^2,...
                Omega_min,Omega_max,optimset('Display','none')); 

            theta(iv_0) = theta_i; 
            Omega(iv_0) = Omega_i;         
            M_g(iv_0)   = Parameter.VSC.NonlinearStateFeedback(Omega_i/Parameter.Turbine.i,theta_i,Parameter); 

        case '2.5' % Determin M_g in Region 2.5, where Omega and theta are fixed
            
            theta_min   = Parameter.CPC.theta_min;
            Omega_rated = Parameter.CPC.Omega_g_rated*Parameter.Turbine.i;
            M_g_min     = 0;
            M_g_max     = Parameter.VSC.M_g_max;
 
            [M_g(iv_0),Omega_dot_Sq(iv_0),exitflag(iv_0)] = ...         
                fminbnd(@(M_g) (OmegaDot(Omega_rated,theta_min,M_g,v_0i,Parameter))^2,...
                M_g_min,M_g_max,optimset('Display','none'));

            theta(iv_0) = theta_min;            
            Omega(iv_0) = Omega_rated;
            
        case '3' % Determine theta in Region 3, where Omega and M_g are fixed   
               
            Omega_i     = Parameter.CPC.Omega_g_rated*Parameter.Turbine.i;
            theta_min   = Parameter.CPC.theta_min;
            theta_max   = max(Parameter.Turbine.SS.theta(:));
            M_g_i       = Parameter.VSC.M_g_rated;

            [theta_i,Omega_dot_Sq(iv_0),exitflag(iv_0)] = ...
                fminbnd(@(theta) (OmegaDot(Omega_i,theta,M_g_i,v_0i,Parameter))^2,...
                theta_min,theta_max,optimset('Display','none'));        
            
            theta(iv_0) = theta_i; 
            Omega(iv_0) = Omega_i;
            M_g(iv_0)   = M_g_i;            
    end
end

%% 5. Calculation of additional variables
lambda   	= Omega*Parameter.Turbine.R./v_0;
c_T      	= interp2(Parameter.Turbine.SS.theta,Parameter.Turbine.SS.lambda,Parameter.Turbine.SS.c_T,theta,lambda,'linear',0);
F_a      	= 1/2*pi*Parameter.Turbine.R^2*c_T*Parameter.General.rho.*v_0.^2;
x_T         = (F_a    )/Parameter.Turbine.k_Te+Parameter.Turbine.x_T0;
P           = M_g*Parameter.Generator.eta_el.*Omega/Parameter.Turbine.i;


%% 6. Plot
figure('Name','Omega')
hold on;grid on;box on;
plot(v_0,radPs2rpm(Omega),'.')
xlabel('v_0 [m/s]')
ylabel('\Omega [rpm]')

figure('Name','theta')
hold on;grid on;box on;
plot(v_0,rad2deg(theta),'.')
xlabel('v_0 [m/s]')
ylabel('\theta [deg]')

figure('Name','M_g')
hold on;grid on;box on;
plot(v_0,M_g,'.')
xlabel('v_0 [m/s]')
ylabel('M_g [Nm]')

figure('Name','x_T')
hold on;grid on;box on;
plot(v_0,x_T,'.')
xlabel('v_0 [m/s]')
ylabel('x_T [m]')

figure('Name','P')
hold on;grid on;box on;
plot(v_0,P,'.')
xlabel('v_0 [m/s]')
ylabel('P [W]')

figure('Name','Torque Controller')
hold on;grid on;box on;
plot(radPs2rpm(Omega),M_g/1e3,'.-')
xlabel('Omega [rpm]')
ylabel('M_g [kNm]')
