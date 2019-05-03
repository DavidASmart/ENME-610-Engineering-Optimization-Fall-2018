
function [C_v] = calc_C_v(D, L_T, C_f, C_p)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the total coefficient of drag of the hull of the AUV (C_v) 
%   given some geometic variables, and the component coefficients of drag
%   due to skin friction and pressure losses (C_f, and C_v)
%   *using the MIT method*
%
%% Inputs:
%   D        outer diameter of the hull (m)
%   L_T      total length of the hull (m)
%   C_f      coefficient of drag of the hull due to skin-friction (n/a - unitless)	
%   C_p      coefficient of drag of the hull due to skin-friction (n/a - unitless)	
%
%% Outputs:
% C_v        total coefficient of drag of the hull (n/a - unitless)
%

%% calculation
C_v = C_f*(1 + (3/2)*(D/L_T)^(3/2) + 7*(D/L_T)^3 + 0.0002*(C_p - 0.6));

end
