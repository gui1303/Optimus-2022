%% StaticCalculationsConfig
% Function: Config for StaticCalculations
%        
% 
%% Usage:
% Adjust CalculationName and run StaticCalculations
%
% Following needs to be declared in the corresponding case:
% - Parameter.Turbine
% Parameter.Turbine.R
% Parameter.Turbine.i
% Parameter.Turbine.SS
% Parameter.Turbine.J
% - Parameter.VSC
% Everything, which is necessary for the VSC 
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
% - check Skywind case
%
%
%% Created:
% David Schlipf on     19-Dec-2014
%
% (c) Universitaet Stuttgart
%

function [v_0,FlagPITorqueControl,Parameter] = StaticCalculationsConfig_Solutions


CalculationName  = 'NREL5MW_FBSWE';       

switch CalculationName
	case {'NREL5MW_FBSWE'}      
        % Case by DS on 02-Dec-2019

  
        %% Default
        Parameter                       	= NREL5MWDefaultParameter_SLOW2DOF;
        Parameter                           = NREL5MWDefaultParameter_FBSWE(Parameter);           

        
        %% NonlinearStateFeedback
        Parameter.VSC.NonlinearStateFeedback    = @(Omega_g,theta,Parameter) min(Parameter.VSC.k*Omega_g^2,Parameter.VSC.M_g_rated); 
        FlagPITorqueControl         	= 1; % 0: only State Feedback, 1: PI controlled in region 1.5 and 2.5
        
        %% Wind speeds
        v_0         = 3.5:.1:30; % [m/s]
        
        %% find v_rated
        v_0_min                         = 0;
        v_0_max                         = 30;        
        Omega                           = Parameter.CPC.Omega_g_rated*Parameter.Turbine.i;
        theta                           = Parameter.CPC.theta_min;  
        M_g                             = Parameter.VSC.M_g_rated;
        [v_rated,Omega_dot_Sq,exitflag] = fminbnd(@(v_0) ...
            (OmegaDot(Omega,theta,M_g,v_0,Parameter))^2,...
            v_0_min,v_0_max,optimset('Display','none'));
        Parameter.VSC.v_rated         	= v_rated;  
        %% find v_1to1d5
        v_0_min                         = 0;
        v_0_max                         = v_rated;        
        Omega                           = Parameter.VSC.Omega_g_1d5*Parameter.Turbine.i;
        theta                           = Parameter.CPC.theta_min;  
        M_g                             = 0;
        [v_1to1d5,Omega_dot_Sq,exitflag]= fminbnd(@(v_0) ...
            (OmegaDot(Omega,theta,M_g,v_0,Parameter))^2,...
            v_0_min,v_0_max,optimset('Display','none'));
        Parameter.VSC.v_1to1d5         	= v_1to1d5;         
        %% find v_1d5to2
        v_0_min                         = v_1to1d5;
        v_0_max                         = v_rated;        
        Omega                           = Parameter.VSC.Omega_g_1d5*Parameter.Turbine.i;
        theta                           = Parameter.CPC.theta_min;  
        M_g                             = Parameter.VSC.Omega_g_1d5^2*Parameter.VSC.k;
        [v_1d5to2,Omega_dot_Sq,exitflag]= fminbnd(@(v_0) ...
            (OmegaDot(Omega,theta,M_g,v_0,Parameter))^2,...
            v_0_min,v_0_max,optimset('Display','none'));
        Parameter.VSC.v_1d5to2         	= v_1d5to2;   
        %% find v_2to2d5
        v_0_min                         = v_1d5to2;
        v_0_max                         = v_rated;        
        Omega                           = Parameter.CPC.Omega_g_rated*Parameter.Turbine.i;
        theta                           = Parameter.CPC.theta_min;  
        M_g                             = Parameter.CPC.Omega_g_rated^2*Parameter.VSC.k;
        [v_2to2d5,Omega_dot_Sq,exitflag]= fminbnd(@(v_0) ...
            (OmegaDot(Omega,theta,M_g,v_0,Parameter))^2,...
            v_0_min,v_0_max,optimset('Display','none'));
        Parameter.VSC.v_2to2d5         	= v_2to2d5;         
    
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
        [v_rated,Omega_dot_Sq,exitflag] = fminbnd(@(v_0) ...
            (OmegaDot(Omega,theta,M_g,v_0,Parameter))^2,...
            v_0_min,v_0_max,optimset('Display','none'));
        Parameter.VSC.v_rated         	= v_rated;        
         

end

end
