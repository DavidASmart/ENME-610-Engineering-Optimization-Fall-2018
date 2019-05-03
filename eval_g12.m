
function [g12] = eval_g12(rho, mu, v, theta, alpha, l, w, d, t, L, f1_g, f1_b, epsilon)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
% evaluates the inequality constraint value g12
% upper limit drag constraint (epsilon constrained / maximizing internal volume)
%
%% Inputs:  
% *parameters*
%   rho        density of water at 25 deg-C (kg/m^3)
%   mu         dynamic viscosity of water at 25 deg-C (N*s/m^2)   
%   v          speed of the AUV (m/s)
%   theta      angle of conic tail section of AUV (deg)
%   alpha      angle of attack of fins (deg)
%   tfins      thickness of hollow fins (m)
% *variables*
%   d          inner diameter of the hull (m)
%   t          thickness of the hull (m)
%   L          length of the cylindrical section of the hull (m)
%   l          length of fins (m)
%   w          width of fins (m)
% *constraint limit*
%   f1_g        "good value" for drag (N)
%   f1_b        "bad value" for drag (N)
%   epsilon     constraint level for drag, as a value between 0 and 1 for
%               the normalized drag
%
%% Outputs:
%   g12        upper bound drag constraint (<= 0)
%

%% prliminary calculations

f1 = eval_f1(rho, mu, v, theta, alpha, l, w, d, t, L);

f1_s = (f1 - f1_g)/(f1_b - f1_g); % normalized

%% inequality constraint
g12 = f1_s/epsilon - 1; % normalized

end