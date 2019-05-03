%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "GridEval"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%       Evaluates the inequality constraints for a full 3D cubic grid of 
%       points spanning the variable bounds. For those points which satisfy 
%       all the constraints, the objective functions are also evaluated.

%
%% Instructions:
%       Just hit "run". It will save the results to GRID_results.mat
%       Plots will be generated showing the feasible domain in the 
%       design space, criterion space, and normalized criterion space.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up

close all
clear
clc

%% get all parameters:
[g, rho, mu, ...
    rho_load, rho_fins, rho_hull,  Sy_hull, ...
    v, depth, T, theta, alpha, tfins, l, w] = set_Params();

%% Variable Bounds:
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generates a set of start points
m = 3;
N = 100;
Nd = N;
Nt = N;
NL = N;
%   m   method of point generation
%         1 - gaussian random
%         2 - halton psudo-random
%         3 - grid
%         4 - sphere
%   N   number of points to generate
%   Nd  number of points in d (for grid only)
%   Nt  number of points in t (for grid only)
%   NL  number of points in L (for grid only)

X0 = genX0(d_L, d_U, t_L, t_U, L_L, L_U, N, Nd, Nt, NL, m);

% plot3(X0(:,1), X0(:,2), X0(:,3),'.')
% axis([d_L, d_U, t_L, t_U, L_L, L_U])
% xlabel('d')
% ylabel('t')
% zlabel('L')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% eval loop
j = 0;
for i = 1:size(X0,1)
    
    d = X0(i,1);
    t = X0(i,2);
    L = X0(i,3);

    [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = eval_gALL(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, theta, alpha, tfins, l, w, ...
    d, t, L, ...
    d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS);

    
    if g1 > 0 || g2  > 0 || g3 > 0 || g4 > 0 || g5 > 0 || g6 > 0 || g7 > 0 || g8 > 0 || g9 > 0 || g10 > 0
        
        % constraints violated
        
    else
        j = j+1;
        
        % state
        X(j,:) = X0(i,:);
        
        % objective
        f1(j,:)	= eval_f1(rho, mu, v, theta, alpha, l, w, d, t, L);  
        f2(j,:)	= eval_f2(theta, d, L);  
    end
    
end

%% good and bad values
f1_b = max(f1);
f1_g = min(f1);

f2_b = min(f2);
f2_g = max(f2);

%% normalize by good and bad values
for i = 1:j
    % scalled values
    f1_s(i,:) = (f1(i) - f1_g)/(f1_b - f1_g);
    f2_s(i,:) = (f2(i) - f2_g)/(f2_b - f2_g);

    % multi-objective Lq-method
    Lq1(i,:) = sum([f1_s(i), f2_s(i)]);
    Lq2(i,:) = sqrt(sumsqr([f1_s(i), f2_s(i)]));
    Lqinf(i,:) = max([f1_s(i), f2_s(i)]);
end

%% save results
save('GRID_results.mat', 'X', 'f1', 'f2', 'f1_s', 'f2_s', 'Lq1', 'Lq2', 'Lqinf')

%%