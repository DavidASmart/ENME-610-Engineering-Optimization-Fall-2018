
function [d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims()
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%% Description:
%   User defines Limit values using this function like a script.
%
%% Inputs:
%   NONE
%
%% Outputs:
%   d_L        lower bound for AUV inner diameter            (m)
%   d_U        upper bound for AUV inner diameter            (m)
%   t_L        lower bound for AUV hull thickness            (m)
%   t_U        upper bound for AUV hull thickness            (m)
%   L_L        lower bound for AUV length                    (m)
%   L_U        upper bound for AUV length                    (m)
%
%   W_lim      upper limit for AUV weight                    (N)
%   FS         factor of safety                              (n/a - unitless)  
%

%% basic variable bounds

% inner-diameter bounds
d_L     = 0.2;      % (m)
d_U     = 0.5;      % (m)

% thickness bounds
t_L     = 0.01;    % (m)
t_U     = 0.05;	   % (m)

% length bounds
L_L     = 1;     % (m)
L_U     = 5;        % (m)

%% constraint bounds

% upper limit weight constraint
W_lim   = 7357.5;   % (N)

% upper limit stress constraint
FS      = 6;        % (n/a - unitless)  

end
