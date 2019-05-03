
function [g, rho, mu, rho_load, rho_fins, rho_hull, Sy_hull, v, depth, T, theta, alpha, tfins, l, w] = set_Params()
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%% Description:
%   User defines Parameter values using this function like a script.
%
%% Inputs:
%   NONE
%
%% Outputs:
%   g          aceleration due to gravity (m/s^2)
%
%   rho        density of water at 25 deg-C (kg/m^3)
%   mu         dynamic viscosity of water at 25 deg-C (N*s/m^2)
%
%   rho_load   density of payload material (kg/m^3)
%   rho_fins   density of fins material (kg/m^3)
%   rho_hull   density of hull material (kg/m^3)
%   Sy_hull    tensile yield strength of hull material (N/m^2)
%
%   v          speed of the AUV (m/s)
%   depth      distance below surface of the water (m)                   
%   T          thrust from the AUV's propeler (N)
%
%   theta      angle of conic tail section of AUV (deg)
%   n          total number of helical ribs  (n/a - unitless) 
%   alpha      angle of attack of fins (deg)
%   tfins      thickness of hollow fins (m)
%   l          length of fins (m)
%   w          width of fins (m)
%

%% gravity
g           = 9.81;      % (m/s^2)

%% properties of water
rho         = 1000;      % (kg/m^3)
mu          = 8.9*10^-3; % (N*s/m^2)

%% propertices of AUV materials

rho_load    = 2070;         % (kg/m^3)
rho_fins    = 1430;         % (kg/m^3) - carbon fiber

% % Al-5058 - (http://asm.matweb.com/)
% rho_hull      = 2660;       % (kg/m^3)
% Sy_hull       = 228*10^6;   % (N/m^2)

% HY-80 Steel - (http://www.matweb.com/)
rho_hull      = 7750;       % (kg/m^3)
Sy_hull       = 552*10^6;   % (N/m^2) 


%% AUV desired specs
v           = 2;        % (m/s)
depth       = 1000;      % (m)

%% AUV hardware
T           = 34.82;    % (N)               <-- Blue Robotics T200 Thruster

%% AUV structural constants
theta       = 30;       % (deg)
alpha       = 5;        % (deg)
tfins       = 0.01;     % (m)
l           = 1;        % (m)
w           = 1;        % (m)

end
