
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "Monotonicity_Analysis_t"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%   Plots the objective functions and inequality constraints as functions
%   of a single variable: t
%
%% Instructions:
%   Just hit "Run". A plot will be generated.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up

close all
clear
clc

% parameters
[g, rho, mu, ...
    rho_load, rho_fins, rho_hull,  Sy_hull, ...
    v, depth, T, theta,alpha, tfins, l, w] = set_Params();

% limits
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();


%% test for a single variable at a time by holding the others constant:

% constant
d = (d_L + d_U)/2;
L = (L_L + L_U)/2;

% varying
N = 100;
t = linspace(t_L, t_U, N);


%% test
for i = 1:N

% Force of Drag (N)
f1(i)	= eval_f1(rho, mu, v, theta, alpha, l, w, d, t(i), L);  

% Internal Volume (m^3)
f2(i)	= eval_f2(theta, d, L); 

% constraints
[g1(i), g2(i), g3(i), g4(i), g5(i), g6(i), g7(i), g8(i), g9(i), g10(i)] = ...
    eval_gALL(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, v, depth, theta, alpha, tfins, l, w,...
    d, t(i), L, ...
    d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS);

end

% rescale / normalize outputs from +/- 0 to 1
f1 = f1/max(abs(f1));
f2 = f2/max(abs(f2));
g3 = g3/max(abs(g3));
g4 = g4/max(abs(g4));
g7 = g7/max(abs(g7));
g8 = g8/max(abs(g8));
g9 = g9/max(abs(g9));
g10 = g10/max(abs(g10));

% plot
figure(1);
hold on
plot(t, f1, 'g', 'LineWidth', 2);
plot(t, f2, 'm', 'LineWidth', 2);
plot(t, g3, 'y', 'LineWidth', 2);
plot(t, g4, 'c', 'LineWidth', 2);
plot(t, g7, 'r', 'LineWidth', 2);
plot(t, g8, 'b', 'LineWidth', 2);
plot(t, g9, 'b', 'LineWidth', 2);
plot(t, g10, 'k', 'LineWidth', 2);
title('Monotonic ?')
xlabel('t (m)');
ylabel('f & g');
legend({'f1','f2','g3','g4','g7','g8','g9','g10'}, 'Location', 'EastOutside');
