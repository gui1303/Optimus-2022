%% StaticCalculationsConfig
% Function: Config for StaticCalculations
%        
% 
%% Usage:
% Adjust CalculationName and run StaticCalculations
%
% Following needs to be declared in the corresponding case:
% - Parameter.Turbine
%   Parameter.Turbine.R
%   Parameter.Turbine.i
%   Parameter.Turbine.SS
%   Parameter.Turbine.J
% - Parameter.VSC
%   Everything, which is necessary for the VSC 
% - v_0: wind speeds to calculate the steady states
% - v_rated
%% Input:
%
% 
%% Output:
% 
%
%% Modified:
%
%
%
%% ToDo:
%
%
%% Created:
% David Schlipf on     19-Dec-2014
%
% (c) Universitaet Stuttgart

function [v_0,FlagPITorqueControl,Parameter] = StaticCalculationsConfig


CalculationName  = 'NREL5MW_FBNREL';       

switch CalculationName
	case {'NREL5MW_FBNREL'}      
        % Case by DS on 24-Nov-2019

  
        %% Default
        Parameter                       	= NREL5MWDefaultParameter_SLOW2DOF;
        Parameter                           = NREL5MWDefaultParameter_FBNREL(Parameter);           

        
        %% NonlinearStateFeedback
     	Parameter.VSC.NonlinearStateFeedback    = @VSControlNREL5MW;        
        FlagPITorqueControl         	= 0; % 0: only State Feedback, 1: PI controlled in region 1.5 and 2.5
        
        %% Wind speeds
        v_0         = 3.5:.1:30; % [m/s]
        
        %% find v_rated
        v_0_min                         = 0;
        v_0_max                         = 30;        
        Omega                           = Parameter.CPC.Omega_g_rated*Parameter.Turbine.i;
        theta                           = Parameter.CPC.theta_min;  
        M_g                             = Parameter.VSC.M_g_rated;
        % Exercise 8.1a: needs adjustments!!!
%         [v_rated,Omega_dot_Sq,exitflag] = fminbnd(@(???) ...
%             ???,...
%             ???,???,optimset('Display','none'));
        Parameter.VSC.v_rated         	= 0;      
        
    case {'NREL5MW_FBSWE'}        
        % Exercise 8.2a: needs adjustments!!!
         

end

end
