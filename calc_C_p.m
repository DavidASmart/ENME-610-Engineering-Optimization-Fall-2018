
function [C_p] = calc_C_p(D, L_T, V)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the coefficient of drag of the hull of the AUV 
%   due to pressure loss (form factor) (C_p) given some geometic variables
%
%% Inputs:
%   D        outer diameter of the hull (m)
%   L_T      total length of the hull (m)
%   V        volume of water displaced by the hull (m^3)
%
%% Outputs:
%   C_p      coefficient of drag of the hull due to skin-friction (n/a - unitless)	
%

%% calculation
C_p = V/(pi*(D/2)^2*L_T);

end
