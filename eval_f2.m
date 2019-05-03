
function [f2] = eval_f2(theta, d, L)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
% evaluates the objective function f2
% Internal Volume of the hull
%
%% Inputs:  
% *variables*
%   d         inner diameter of the hull (m)
%   L         length of the hull (m)
%
%% Outputs:
%   f2        Internal Volume of Cylindrical Section of the hull (m^3)
%

%% objective function
f2 = calc_V(theta, d, L);

end
