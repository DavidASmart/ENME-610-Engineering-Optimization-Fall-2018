
function [s_max] = calc_s_max(s_t, s_h, s_a, s_HH, s_AA)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the Von-Mises combined stresses in the AUV sections
%   ... and the maximum throughout the entire AUV hull
%   based on the component stresses
%
%% Inputs:
%   s_t      max tangential stresses in a sphere (N)
%   s_h      max hoop-stress in a cylinder (N)
%   s_a      max axial stress in a cylinder (N)
%   s_HH     max hoop stress in a cone (N)
%   s_AA     max axial stress in a cone (N)
%
%% Outputs:
%   s_sphere       Von-Mises Stress in the spherical front end of the hull (N)
%   s_cylinder     Von-Mises Stress in the cylindrical section of the hull (N)
%   s_cone         Von-Mises Stress in the conical tail section of the hull (N)
%   s_max          maximum Von-Mises stress in the entire AUV hull (N)
%

%% calculation

% front hemispherical section
s_sphere = calc_s_sphere(s_t);    
        
% middle cylindrical section
s_cylinder = calc_s_cylinder(s_h, s_a);

% conic tail section
s_cone = calc_s_cone(s_HH, s_AA);

% max of all sections
s_max = max([s_sphere, s_cylinder, s_cone]);

end
