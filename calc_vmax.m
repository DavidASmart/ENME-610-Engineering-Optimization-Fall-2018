
function [v] = calc_vmax(g, rho, mu, D, L_T, V, S, A, T, C_D_fins, W)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Calculates the Terminal Velocity of AUV, 
%   This function was created just to check that our calculations were 
%   correct, and confirmed them. It showed that the AUV could reach or 
%   exceed the set parameter value for velocity (v) when the drag 
%   calculated at v was less than the available thrust T.
%   This function could potetentially be used as another objective in a 
%   future iteration of this project.
%
%   Also, currently set up to make a plot for drag force, net forward force, 
%   and speed vs time.
%
%% Inputs:
%   g          aceleration due to gravity               (m/s^2)
%   rho        density of water at 25 deg-C             (kg/m^3)
%   mu         dynamic viscosity of water at 25 deg-C   (N*s/m^2)
%   D          outer diameter of the hull               (m)
%   L_T        total length of the hull                 (m)
%   V          volume of water displaced by the hull    (m^3)
%   S          surface area of the hull                 (m^2)
%   A          cross-section area of fins               (m^2)
%   T          thrust from the AUV's propeler           (N)
%   C_D_fins   Drag Coefficient for fins                (n/a - unitless)
%   W          weight of the empty AUV                  (N)
%
%% Outputs:
%   v          terminal velocity                        (m/s)
%

a = 1; % just to get into the loop

% initialize
F_D = 0; v = 0;
t = 0; dt = 0.01;

% stopping criteria
tol = 1*10^-6;

figure; hold on

while a > tol
    % calc thrust - drag
    F = T - F_D;

    % calc acceleration
    a = F/(W/g);

    % calc velocity
    v = v + a*dt;

    % reynold's number
    Rn = calc_Rn(rho, mu, v, L_T);

    % Drag Coefficient for hull (n/a - unitless)
    [C_v] = calc_HydroCoeff(D, L_T, V, Rn);
    
    % calc drag
    F_D_hull = C_v*(rho/2)*S*v^2;
    F_D_fin = C_D_fins*(rho/2)*A*v^2;
    F_D = F_D_hull + F_D_fin;
    
    t = t+dt;
    
    % plot
    plot(t, F_D,'r.')
    plot(t, F,'g.')
    plot(t, v,'b.')
end

pause(0.1)

end