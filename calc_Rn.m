
function [Rn] = calc_Rn(rho, mu, v, L_T)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the Reynolds number (Rn) given the AUV's total length and  velocity
%
%% Inputs:
%   rho        density of water at 25 deg-C (kg/m^3)
%   mu         dynamic viscosity of water at 25 deg-C (N*s/m^2)
%   v          speed of the AUV (m/s)
%   L_T        total length of the hull (m)
%
%% Outputs:
%   Rn         Reynold's Number (n/a - unitless)
%

%% calculation
Rn = (rho*v*L_T)/mu;

end
