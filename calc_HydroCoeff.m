
function [C_v] = calc_HydroCoeff(D, L_T, V, Rn)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the total coefficient of drag of the hull of the AUV (C_v)
%   based on geometric quantities and Reynolds Number (Rn)
%
%% Inputs:
%   D        outer diameter of the hull (m)
%   L_T      total length of the hull (m)
%   V        volume of water displaced by the hull (m^3) 	
%   Rn       Reynold's Number (n/a - unitless)
%
%% Outputs:
%   C_v      total coefficient of drag of the hull (n/a - unitless)
%

%% calculation
C_f = calc_C_f(Rn);                 % coefficient of drag due to skin-friction (n/a - unitless)	
C_p = calc_C_p(D, L_T, V);          % coefficient of drag due to pressure loss/form factor (n/a - unitless)	
C_v = calc_C_v(D, L_T, C_f, C_p);   % total coefficient of drag (n/a - unitless)	

end
