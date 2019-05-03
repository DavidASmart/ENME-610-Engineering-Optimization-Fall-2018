
function [s_sphere] = calc_s_sphere(s_t)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the Von-Mises combined stresses in the spherical front section of the hull
%   based on the component stresses
%
%% Inputs:
%   s_t      max tangential stresses in a sphere (N)
%
%% Outputs:
%   s_sphere       Von-Mises Stress in the spherical front end of the hull (N)
%

%% calculation
s_sphere = calc_vonMises(s_t, s_t, 0);

end
