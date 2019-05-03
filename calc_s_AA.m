
function [s_AA] = calc_s_AA(D, d, theta, P)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the axial stress in the conical tail section of the AUV
%   *Using the thiin-wall aproximate pressure vessel equation*
%
%% Inputs:
%   D        outer diameter of the hull (m)
%   d        inner diameter of the hull (m)
%   theta    angle of the conic tail section of the hull (deg)
%   P        hydrostatic pressure (N)
%
%% Outputs:
%   s_AA      max axial stress in a cone (N)
% 

%% calculation
s_AA = P*(1/2)*D*sind(theta)*tand(theta)/(D - d);

end
