
function [L_T] = calc_L_T(theta, D, L)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the total length of the AUV (L_T) given the variables D, L, and theta
%
%% Inputs:
%   theta      angle of the conic tail section of the hull (deg)
%
%   D          outer diameter of the hull (m)
%   L          length of the cylindrical section of the hull (m)
%
%% Outputs:
%   L_T        total length of the hull	(m)
%

%% calculation
L_T = D/2 + L + (D/2)*(1/tand(theta));

end
