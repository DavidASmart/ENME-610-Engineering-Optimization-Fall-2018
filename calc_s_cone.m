
function [s_cone] = calc_s_cone(s_HH, s_AA)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the Von-Mises combined stresses in the conical tail section of the AUV hull 
%   based on the component stresses
%
%% Inputs:
%   s_HH      max hoop stress in a cone (N)
%   s_AA      max axial stress in a cone (N)
%
%% Outputs:
%   s_cone  Von-Mises Stress in the conical tail section of the hull (N)
%

%% calculation
s_cone = calc_vonMises(s_HH, s_AA, 0);

end
