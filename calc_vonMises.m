
function [s_vonMises] = calc_vonMises(s_1, s_2, s_3)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the Von-Mises combined stresses
%   based on the principle stresses
%
%% Inputs:
%   s_1      1st principle stress (N/m^2)
%   s_2      2nd principle stress (N/m^2)
%   s_3      3rd principle stress (N/m^2)
%
%% Outputs:
%   s_vonMises     Von-Mises Max Stress(N/m^2)
%

%% calculation
s_vonMises = (1/sqrt(2))*sqrt((s_1 - s_2)^2 + (s_1 - s_3)^2 + (s_3 - s_2)^2);

end
