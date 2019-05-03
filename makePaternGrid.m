
function [X0] = makePaternGrid(d_L, d_U, t_L, t_U, L_L, L_U, Nd, Nt, NL)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
% makes a 3D cubic grid of points cover the full range of varibale values
% with resolution of points dictated by Nd, Nt, and NL
%
%% Inputs:
%   d_L        lower bound for AUV inner diameter            (m)
%   d_U        upper bound for AUV inner diameter            (m)
%   t_L        lower bound for AUV hull thickness            (m)
%   t_U        upper bound for AUV hull thickness            (m)
%   L_L        lower bound for AUV length                    (m)
%   L_U        upper bound for AUV length                    (m)
%   Nd         number of points in d (for grid only)
%   Nt         number of points in t (for grid only)
%   NL         number of points in L (for grid only)
%
%% Outputs:
%   X0          set of start points
%

%% make grid
dd = linspace(d_L, d_U, Nd);
tt = linspace(t_L, t_U, Nt);
LL = linspace(L_L, L_U, NL);
[ddd, ttt, LLL] = ndgrid(dd, tt, LL);

%% reformat
X0 = [ddd(:)';ttt(:)';LLL(:)']';

end
