
function [V_fins2] = calc_V_fins2(tfins, l, w)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the mass-volume of the hollow fins
%
%% Inputs:
%   tfins    thickness of fins material (m)
%
%   l        length of fins (m)
%   w        width of fins (m)
%
%% Outputs:
%   V_fins2   volume of water displaced by the fins (m^3)
%

%% calculation
cross_section = (pi/4)*((0.12*l)*l - (0.12-2*tfins)*l*(l-2*tfins));
end_cap = (pi/4)*l*(0.12*l)*tfins;
V_fins2 = 2*cross_section*(w-tfins) + 2*end_cap; % http://structx.com/Shape_Formulas_035.html

end
