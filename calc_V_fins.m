
function [V_fins] = calc_V_fins(l, w)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the volume of water displaced by the fins
%
%% Inputs:
%   l        length of fins (m)
%   w        width of fins (m)
%
%% Outputs:
%   V_fins   volume of water displaced by the fins (m^3)
%

%% calculation
cross_section = (pi/4)*l*(0.12*l); % approx. as ellispe (NACA 0012 means 0.12 height to length ration)
V_fins = 2*(cross_section*w);

end
