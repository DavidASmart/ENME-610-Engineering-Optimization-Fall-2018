%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "CheckFeas"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Checks optima for feasibility
%
%% Instructions:
%   Fill in the variable x with the optima point, then hit "Run"
%   It will display the values for each of the inequality constraints.
%   It will display Y/N depending on if the point is feasible.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% clean up
close all
clear
clc

%% point to be evaluated

x = [x1, x2, x3]; % to be filled in by user***

d = x(1);
t = x(2);
L = x(3);

%% get all parameters:
[g, rho, mu, ...
    rho_load, rho_fins, rho_hull,  Sy_hull, ...
    v, depth, T, theta, alpha, tfins, l, w] = set_Params();

%% Variable Bounds:
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

%% check feasibility
[g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = ...
    eval_gALL(g, rho, rho_load, rho_fins, rho_hull, Sy_hull, v, depth, theta, alpha, tfins, l, w, d, t, L, d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS)

if g1 > 0 || g2  > 0 || g3 > 0 || g4 > 0 || g5 > 0 || g6 > 0 || g7 > 0 || g8 > 0 || g9 > 0 || g10 > 0
    fprintf('\n\n N \n\n')
else
    fprintf('\n\n Y \n\n')
end

