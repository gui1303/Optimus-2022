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

%Gain scheduling with 0.00150 minimum and 0.0040 max
RefGain  = [0.0010 0.0060];
RefTheta = [min(theta) max(theta)]; 
kp       = interp1(RefTheta,RefGain,theta);
Ti       = 10;
ki       = kp./Ti;

fprintf('theta      = [%s];\n',sprintf('%f ',theta));
fprintf('kp         = [%s];\n',sprintf('%f ',-kp));
fprintf('ki         = [%s];\n',sprintf('%f ',-ki));