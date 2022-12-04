% -----------------------------
% Script: Finds rated wind speed
% Exercise 08 of Master Course 
% "Controller Design for Wind Turbines and Wind Farms"
% ------------
% Task:
% - adjust script (fminbnd,fminunc)
% ------------
% History:
% v01:	David Schlipf on 22-Nov-2021
% ----------------------------------
clear all;clc;close all;
Parameter                               = NREL5MWDefaultParameter_SLOW2DOF;
Parameter                               = NREL5MWDefaultParameter_FBNREL(Parameter);   
v_0                                     = 5:.1:30; % [m/s]
v_0_min                                 = 0;
v_0_max                                 = 30;        
Omega                                   = Parameter.CPC.Omega_g_rated*Parameter.Turbine.i;
theta                                   = Parameter.CPC.theta_min;  
M_g                                     = Parameter.VSC.M_g_rated;
        
%% Brute Force Optimization       
for  iv_0=1:length(v_0)        
    Residual(iv_0)=OmegaDot(Omega,theta,M_g,v_0(iv_0),Parameter); 
end   

fun = @(x)(OmegaDot(Omega,theta,M_g,x,Parameter))^2; 
minimobnd = fminbnd(fun,v_0_min,v_0_max);
minimounc = fminunc(fun,10);

figure
hold on
plot(minimobnd,0,'o')
plot(v_0,Residual*60/2/pi)
plot([v_0(1) v_0(end)],[0,0])
xlabel('wind speed [m/s]')
ylabel('rotor acceleration [rpm/s]')


%% Optimization using fminbnd
% [v_rated,Omega_dot_Sq,exitflag] = fminbnd(...,optimset('Display','iter'));

%% Optimization using fminunc
% [v_rated,Omega_dot_Sq,exitflag] = fminunc(...,optimset('Display','iter'));