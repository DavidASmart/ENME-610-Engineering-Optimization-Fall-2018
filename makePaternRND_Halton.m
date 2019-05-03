
function [X0] = makePaternRND_Halton(d_L, d_U, t_L, t_U, L_L, L_U, N)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
% samples "N" psudo-random points from a halton set 
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
h = haltonset(3);

%% scale to span variable bounds
LB = [d_L, t_L, L_L]; 
UB = [d_U, t_U, L_U]; 
X0 = LB + h(1:N,:).*(UB - LB);

end
