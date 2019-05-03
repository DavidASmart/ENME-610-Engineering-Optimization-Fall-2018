
function [V_hull] = calc_V_hull(theta, D, d, L)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the mass-volume of the AUV (V_hull) given d, D, L, and theta
%
%% Inputs:
%   theta      angle of the conic tail section of the hull (deg)
%
%   D          outer diameter of the hull (m)
%   d          inner diameter of the hull (m)
%   L          length of the cylindrical section of the hull (m)
%
%% Outputs:
%   V_hull     volume of the hull material (m^3)
%

%% calculation
V_hemisphere = (1/12)*pi*(D^3 - d^3);           % volume of front hemisphere
V_cylinder = (1/4)*pi*(D^2 - d^2)*L;            % volume of middle cylinder
V_cone = (1/24)*pi*(D^3 - d^3)*(1/tand(theta)); % volume of connic tail
V_hull = V_hemisphere + V_cylinder + V_cone;    % sum total of all

end
