
function [s_a] = calc_s_a(D, d, P)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the axial stress in the cylindrical middle section of the AUV
%   *Using the thick-wall full pressure vessel equation*
%
%% Inputs:
%   D       outer diameter of the hull (m)
%   d       inner diameter of the hull (m)
%   P       hydrostatic pressure (N)
%
%% Outputs:
%   s_a    max axial stress in a cylinder (N)
% 

%% calculation
s_a = P*d^2/(D^2 - d^2);

end
