
function [s_t, s_h, s_a, s_HH, s_AA] = calc_PV_Stresses(theta, d, D, P)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the pressure-vessel component stresses in the AUV hull based
%   on the basic geometric properties, and of course hydrostatic pressure
%   *Using the thick-wall full pressure vessel equation*
%
%% Inputs:
%   theta	 angle of the conic tail section of the hull (deg)
%
%   d        inner diameter of the hull (m)
%   D        outer diameter of the hull (m)
%   P        hydrostatic pressure (N)
%
%% Outputs:
%   s_t      max tangential stresses in a sphere (N)
%   s_h      max hoop-stress in a cylinder (N)
%   s_a      max axial stress in a cylinder (N)
%   s_HH     max hoop stress in a cone (N)
%   s_AA     max axial stress in a cone (N)

%% calculation

% front hemispherical section
s_t = calc_s_t(D, d, P);

% middle cylindrical section
s_h = calc_s_h(D, d, P);
s_a = calc_s_a(D, d, P);

% conic tail section
s_HH = calc_s_HH(D, d, theta, P);
s_AA = calc_s_AA(D, d, theta, P);

end
