
function [C_f] = calc_C_f(Rn)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the coefficient of drag of the hull of the AUV 
%   due to skin-friction (C_f) given the Reynold's Number
%
%% Inputs:
%   Rn        Reynold's Number (n/a - unitless)
%
%% Outputs:
%   C_f       coefficient of drag of the hull due to skin-friction (n/a - unitless)	
%

%% calculation
C_f = 0.0075/(log10(Rn)-2)^2;

end
