
function [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = ...
    eval_gALL(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, v, depth, theta, alpha, tfins, l, w, ...
    d, t, L, ...
    d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS)
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
% evaluates all the inequality constraints values
%
%% Inputs:  
% *parameters*
%   g          aceleration due to gravity                    (m/s^2)
%   rho        density of water at 25 deg-C                  (kg/m^3)
%   rho_load   density of payload material                   (kg/m^3)
%   rho_fins   density of fins material                      (kg/m^3)
%   rho_hull   density of hull material                      (kg/m^3)
%   Sy_hull    tensile yield strength of hull material       (N/m^2)
%   v          speed of the AUV                              (m/s)
%   depth      distance below surface of the water           (m)         
%   theta      angle of conic tail section of AUV            (deg)
%   alpha      angle of attack of fins                       (deg)
%   tfins      thickness of hollow fins                      (m)
%   l          length of fins                                (m)
%   w          width of fins                                 (m)
% *variables*
%   d          inner diameter of the hull                    (m)
%   t          thickness of the hull                         (m)
%   L          length of the cylindrical section of the hull (m)
% *constraint limit*
%   d_L        lower bound for AUV inner diameter            (m)
%   d_U        upper bound for AUV inner diameter            (m)
%   t_L        lower bound for AUV hull thickness            (m)
%   t_U        upper bound for AUV hull thickness            (m)
%   L_L        lower bound for AUV length                    (m)
%   L_U        upper bound for AUV length                    (m)
%   W_lim      upper limit for AUV weight                    (N)
%   FS         factor of safety                              (n/a - unitless)    
%
%% Outputs:
%   g1         inner-diameter lower bound constraint
%   g2         inner-diameter upper bound constraint
%   g3         thickness lower bound constraint
%   g4         thickness upper bound constraint
%   g5         length lower bound constraint
%   g6         length upper bound constraint
%   g7        upper bound weight constraint
%   g8        lower bound bouyancy constraint
%   g9        upper bound bouyancy constraint
%   g10        upper bound stress constraint
%

%% prliminary calculations

% outer diameter
D = d + 2*t;

% volume displaced
V = calc_V(theta, D, L);

% volume displaced by fins
V_fins = calc_V_fins(l, w);

% cross-section area of fins
A = 2*w*l;

% bouyancy
F_B = rho*g*(V + V_fins);

% Lift Coefficient for fins
C_L_fins = 0.1089*alpha + 0.0194;

% force of lift from fins
F_L = C_L_fins*(rho/2)*A*v^2;

% mass-volume of hull
V_hull = calc_V_hull(theta, D, d, L);

% mass-volume of fins
V_fins2 = calc_V_fins2(tfins, l, w);

% weight
W_hull = g*rho_hull*V_hull;
W_fins = g*rho_fins*V_fins2;
W = W_hull + W_fins;

% Internal Volume of Cylindrical Section of the hull (m^3)
V_i = calc_V(theta, d, L);

% resulting payload weight
W_load = g*rho_load*V_i;

% stress
P = rho*g*depth; % hydrostatic pressure
[s_t, s_h, s_a, s_HH, s_AA] = calc_PV_Stresses(theta, d, D, P);
s_max = calc_s_max(s_t, s_h, s_a, s_HH, s_AA);

%% inequality constraints

% inner-diameter bounds
% g1 = d_L - d;
% g2 = d - d_U;
g1 = 1 - d/d_L; % normalized
g2 = d/d_U - 1; % normalized

% thickness bounds
% g3 = t_L - t;
% g4 = t - t_U;
g3 = 1 - t/t_L; % normalized
g4 = t/t_U - 1; % normalized

% length bounds
% g5 = L_L - L;
% g6 = L - L_U;
g5 = 1 - L/L_L; % normalized
g6 = L/L_U - 1; % normalized

% upper limit weight constraint
% g7 = W + W_load - W_lim;
g7 = (W + W_load)/W_lim - 1; % normalized

% net bouyancy bounds
% g8 = (W + W_load) - F_B - F_L;
% g9 = F_B - (W + W_load) - F_L;
g8 = (W + W_load - F_B)/F_L - 1; % normalized
g9 = (F_B - W - W_load)/F_L - 1; % normalized

% upper limit stress constraint
% g10 = s_max - (1/FS)*Sy_hull;
g10 = s_max/((1/FS)*Sy_hull) - 1; % normalized

end
