
function [X0] = makePaternRND(d_L, d_U, t_L, t_U, L_L, L_U, N)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
% samples "N" random points from a gaussian distribution 
% covering the full range of varibale values
%
%% Inputs:
%   d_L        lower bound for AUV inner diameter            (m)
%   d_U        upper bound for AUV inner diameter            (m)
%   t_L        lower bound for AUV hull thickness            (m)
%   t_U        upper bound for AUV hull thickness            (m)
%   L_L        lower bound for AUV length                    (m)
%   L_U        upper bound for AUV length                    (m)
%   N          number of points to generate
%
%% Outputs:
%   X0          set of start points
%

%% generate points
rng(1)
X00 = rand(N,3);

%% scale points to span the variable bounds
LB = [d_L, t_L, L_L]; 
UB = [d_U, t_U, L_U]; 

X0 = LB + X00.*(UB - LB);

end
