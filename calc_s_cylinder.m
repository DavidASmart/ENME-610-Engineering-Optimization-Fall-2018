
function [s_cylinder] = calc_s_cylinder(s_h, s_a)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the Von-Mises combined stresses in the cylindrical middle section of the AUV hull 
%   based on the component stresses
%
%% Inputs:
%   s_h      max hoop-stress in a cylinder (N)
%   s_a      max axial stress in a cylinder (N)
%
%% Outputs:
%   s_cylinder     Von-Mises Stress in the cylindrical section of the hull (N)
%

%% calculation
s_cylinder = calc_vonMises(s_h, s_a, 0);

end
