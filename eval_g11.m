
function [g11] = eval_g11(d, L, f2_g, f2_b, epsilon)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
% evaluates the inequality constraint value g11
% lower limit interal volume constraint (epsilon constrained / minimizing drag)
%
%% Inputs:  
% *variables*
%   d          inner diameter of the hull (m)
%   L          length of the cylindrical section of the hull (m)
%
% *constraint limit*
%   f2_g       "good value" for internal volume (m^3)
%   f2_b       "bad value" for internal volume (M^3)
%   epsilon    constraint level for volume, as a value between 0 and 1 for
%              the normalized volume
%
%% Outputs:
%   g11        lower bound interal volume constraint (<= 0)
%

%% preliminary
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, theta, ~, ~, ~, ~] = set_Params();

%% objective function
f2 = eval_f2(theta, d, L);

f2_s = (f2 - f2_g)/(f2_b - f2_g); % normalized

%% inequality constraint
g11 = f2_s/epsilon - 1; % normalized

end
