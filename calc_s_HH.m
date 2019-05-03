
function [s_H] = calc_s_HH(theta, D, d, P)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the hoop stress in the conical tail section of the AUV
%   *Using the thin-wall aproximate pressure vessel equation*
%
%% Inputs:
%   theta    angle of the conic tail section of the hull (deg)
%
%   D        outer diameter of the hull (m)
%   d        inner diameter of the hull (m)
%   P        hydrostatic pressure (N)
%
%% Outputs:
%   s_H      max hoop stress in a cone (N)
% 

%% calculation
s_H = P*D*sind(theta)*tand(theta)/(D - d);

end
