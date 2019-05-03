
function [V] = calc_V(theta, D, L)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the volume of water displaced by the AUV (V) given the variables D, L, and theta
%
%% Inputs:
%   theta      angle of the conic tail section of the hull (deg)
%
%   D          outer diameter of the hull (m)
%   L          length of the cylindrical section of the hull (m)
%
%% Outputs:
%   V          volume of water displaced by the hull (m^3)
%

%% calculation
V = (1/12)*pi*D^3 + (1/4)*pi*D^2*L + (1/24)*pi*D^3*(1/tand(theta));

end
