
function [s_t] = calc_s_t(D, d, P)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the tangential stress in the hemispherical front section of the AUV
%   *Using the thick-wall full pressure vessel equation*
%
%% Inputs:
%   D        outer diameter of the hull (m)
%   d        inner diameter of the hull (m)
%   P        hydrostatic pressure (N)
%
%% Outputs:
%   s_t      max tangential stresses in a sphere (N) 
% 

%% calculation
s_t = P*D^2*2/(D^2-d^2);

end
