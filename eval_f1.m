
function [f1] = eval_f1(rho, mu, v, theta, alpha, l, w, d, t, L)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
% evaluates the objective function f1
% Force of Drag on the AUV
%
%% Inputs:  
% *parameters*
%   rho        density of water at 25 deg-C (kg/m^3)
%   mu         dynamic viscosity of water at 25 deg-C (N*s/m^2)   
%   v          speed of the AUV (m/s)
%   theta      angle of conic tail section of AUV (deg)
%   alpha      angle of attack of fins (deg)
%   l          length of fins (m)
%   w          width of fins (m)
% *variables*
%   d          inner diameter of the hull (m)
%   t          thickness of the hull (m)
%   L          length of the cylindrical section of the hull (m)
%
%% Outputs:
%   f1        Force of Drag on the AUV (N)
%

%% prliminary calculations

% outer diameter
D = d + 2*t;

% total length
L_T = calc_L_T(theta, D, L);

% volume displaced
V = calc_V(theta, D, L);

% surface area of hull
S = calc_S(theta, D, L);

% cross-section area of fins
A = 2*w*l;

% reynold's number
Rn = calc_Rn(rho, mu, v, L_T);

% Drag and Lift Coefficients (n/a - unitless)
[C_v] = calc_HydroCoeff(D, L_T, V, Rn);

% Drag and Lift Coefficients for fins (n/a - unitless)
C_D_fins = 0.00012*alpha^2 + 0.0004*alpha + 0.0006;

% force of drag on the hull (N)
F_D = C_v*(rho/2)*S*v^2;

% force of drag on the fins (N)
F_D_fins = C_D_fins*(rho/2)*A*v^2;

% total force of drag on the AUV (N)
F_D = F_D + F_D_fins;

%% objective function
f1 = F_D;

end