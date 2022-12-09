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
            Omega(iv_0) = NaN;
            M_g(iv_0)   = NaN;
            theta(iv_0) = NaN;
        
                   
        case '3' % Determin theta in Region 3, where Omega and M_g are fixed   
            % Exercise 8.1b: needs adjustments!!!   
            Omega(iv_0) = NaN;
            M_g(iv_0)   = NaN;
            theta(iv_0) = NaN;           
    end
end

%% 5. Calculation of additional variables
% Exercise 8.1b: needs adjustments!!!
x_T         = NaN;
P           = NaN;

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
