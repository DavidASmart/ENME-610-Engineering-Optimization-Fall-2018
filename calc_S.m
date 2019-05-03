
function [S] = calc_S(theta, D, L)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the surface area of the AUV (S) given the variables D, L, and theta
%
%% Inputs:
%   theta      angle of the conic tail section of the hull (deg)
%
%   D          outer diameter of the hull (m)
%   L          length of the cylindrical section of the hull (m)
%
%% Outputs:
%   S          surface area of the hull (m^2)

%% calculation
S = (1/2)*pi*D^2 + pi*D*L + (1/4)*pi*D^2*abs(1/sind(theta));

end
