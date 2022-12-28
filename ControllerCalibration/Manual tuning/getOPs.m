%pitch controller gain scheduling
clear all;clc;close all
OPs     = [12 14 18 24];
nOP     = length(OPs);
theta   = NaN(nOP,1);
SteadyStates = load('SteadyStatesOptimusFAST.mat');
for iOP=1:nOP  
    OP         = OPs(iOP);
    index      = find(abs(SteadyStates.v_0 - OP) < 0.001);
    theta(iOP) = SteadyStates.theta(index);
end

theta(end + 1) = deg2rad(90);

%Gain scheduling 
RefGain  = [0.001400 0.005];
RefTheta = [min(theta) max(theta)]; 
kp       = interp1(RefTheta,RefGain,theta);
Ti       = 10;
ki       = kp./Ti;

fprintf('theta      = [%s];\n',sprintf('%f ',theta));
fprintf('kp         = [%s];\n',sprintf('%f ',-kp));
fprintf('ki         = [%s];\n',sprintf('%f ',-ki));