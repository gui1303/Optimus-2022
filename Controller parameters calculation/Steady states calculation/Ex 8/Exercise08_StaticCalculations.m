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
[v_0,FlagPITorqueControl,Parameter]    = StaticCalculationsConfig;

OmegaRated  = 0.8221;
OmegaMax    = Parameter.Turbine.i*Parameter.VSC.Omega_g_2_5To3;
OmegaMin    = Parameter.Turbine.i*Parameter.VSC.Omega_g_1To1_5;
ThetaMax    = max(Parameter.Turbine.SS.theta);
ThetaMin    = min(Parameter.Turbine.SS.theta);

%% 3. Allocation
Omega               = zeros(1,length(v_0));
theta               = zeros(1,length(v_0));
M_g                 = zeros(1,length(v_0));
Omega_dot_Sq        = zeros(1,length(v_0));
exitflag            = zeros(1,length(v_0));

%% 4. Loop over wind speeds to determine omega, theta, M_g
for iv_0=1:length(v_0)
    v_0i        	= v_0(iv_0);
    
    %% 4.1 Determin Region
    if FlagPITorqueControl
            % Exercise 8.2a: needs adjustments!!!
      
    else % no PI torque control
        if      v_0i < Parameter.VSC.v_rated
            Region = 'StateFeedback';
        else
            Region = '3';
        end
    end

    %% 4.2 Determin Static Values
    switch Region %
        case {'2','StateFeedback'} % Determin Omega and M_g in Region 2 (or 1-2.5 for state feedback), where theta is fixed 
            % Exercise 8.1b: needs adjustments!!!
            theta(iv_0) = Parameter.CPC.theta_min;                     % theta min in rad
            fun         = @(x)(OmegaDot(x,theta(iv_0),Parameter.VSC.k*(x/Parameter.Turbine.i)^2,v_0i,Parameter))^2;
            Omega(iv_0) = (fminbnd(fun,OmegaMin,OmegaMax));
            M_g(iv_0)   = Parameter.VSC.NonlinearStateFeedback(Omega(iv_0)/Parameter.Turbine.i,theta(iv_0),Parameter);    
        case '3' % Determin theta in Region 3, where Omega and M_g are fixed   
            % Exercise 8.1b: needs adjustments!!!  
            Omega(iv_0) = OmegaRated;
            M_g(iv_0)   = Parameter.VSC.M_g_rated;
            fun         = @(x)(OmegaDot(Omega(iv_0),x,M_g(iv_0),v_0i,Parameter))^2;
            theta(iv_0) = fminbnd(fun,ThetaMin,ThetaMax);         
    end
end

%% 5. Calculation of additional variables
% Exercise 8.1b: needs adjustments!!!
x_T         = NaN;
P = NaN;
Fa = NaN;
for iv_0=1:length(v_0)
    v_0i        	= v_0(iv_0);
    Fa(iv_0)  = AerodynamicThrust(Omega(iv_0),theta(iv_0),v_0i,Parameter);
    x_T(iv_0) = Fa(iv_0)/Parameter.Turbine.k_Te + Parameter.Turbine.x_T0;
    P(iv_0)   = M_g(iv_0)*(Omega(iv_0)/Parameter.Turbine.i);
end

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

function F_a = AerodynamicThrust(Omega,theta,v,Parameter)
    lambda      = Omega*Parameter.Turbine.R/v;
    c_T         = interp2(Parameter.Turbine.SS.theta,Parameter.Turbine.SS.lambda,Parameter.Turbine.SS.c_T,theta,lambda,'linear',0);
    F_a         = (1/2*pi*Parameter.Turbine.R^3*Parameter.General.rho*c_T/lambda*v^2);
end
