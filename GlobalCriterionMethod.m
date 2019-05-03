
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "GlobalCriterionMethod"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%       A custom Lq method multi-objective optimizer (utilizing fmincon for
%       constrained single-objective optimization).
%% Instructions:
%       Choose a value for q, and then hit "Run". It will print out the
%       chosen q the multi-objective solution.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up

close all
clear
clc

% ** q = 1, 2, inf **
q = 2;

% limits
[~, rho, mu, rho_load, rho_fins, rho_hull, Sy_hull, v, depth, T, theta, alpha, tfins, l, w] = set_Params();
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

% variable bounds
LB = [d_L, t_L, L_L]; % Variable Lower Bounds
UB = [d_U, t_U, L_U]; % Variable Upper Bounds

% initial guess
d0 = d_L + (d_U - d_L)*(1);
t0 = t_L + (t_U - t_L)*(1);
L0 = L_L + (L_U - L_L)*(1);
X0 = [d0, t0, L0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimization / minimization

% currently these options are pretty much default, but I put this here for
% potential uses in the future
options = optimoptions(@fmincon, ...
    'ConstraintTolerance',      1e-6, ...
    'MaxFunctionEvaluations',   2000, ...
    'MaxIterations',            1500, ...
    'StepTolerance',            1e-10);

% single objective optimization to find the ideal point
[X_opt1, f1_0] = fmincon(@(X)myfun1(X), X0, [], [], [], [], LB, UB, @mycon, options);
[X_opt2, f2_0] = fmincon(@(X)myfun2(X), X0, [], [], [], [], LB, UB, @mycon, options);

% multi-objective optimization according to the Lq method
[X_opt, f_opt] = fmincon(@(X)myfunLq(X, q, f1_0, f2_0), X0, [], [], [], [], LB, UB, @mycon, options);

%% minimized objective values

d = X_opt(1);
t = X_opt(2);
L = X_opt(3);

% Force of Drag (N)
f1	= eval_f1(rho, mu, v, theta, alpha, l, w, d, t, L);  

% Internal Volume (m^3)
f2	= eval_f2(theta, d, L); 

% good and bad values
f1_g    = 23.6389;
f1_b    = 25.1932;
f2_b    = 0.0353;
f2_g    = 0.1242;

% scalled values
f1_s = (f1 - f1_g)/(f1_b - f1_g);
f2_s = (f2 - f2_g)/(f2_b - f2_g);

%% display results
fprintf('\n\n Optimization Method: ')
fprintf('\n Global Criterion (Lq) Method [with fmincon] ')

fprintf('\n\n Optimizer Parameters: ')
fprintf('\n q = %3.3f ', q);

fprintf('\n\n Multi-Objective Optimum: ')
fprintf('\n D = %3.3f m ', d);
fprintf('\n t = %3.3f m ', t);
fprintf('\n L = %3.3f m ', L);

fprintf('\n\n Optimized Objective:')
fprintf('\n f1 = %3.3f N', f1)
fprintf('\n f2 = %3.3f m^3', f2)

fprintf('\n\n Scaled Values: ')
fprintf('\n f1s = %3.3f', f1_s)
fprintf('\n f2s = %3.3f', f2_s)
fprintf('\n\n')

%% save results
save(strcat('Lq',num2str(q),'_results.mat'), 'X_opt', 'f1', 'f2', 'f1_s', 'f2_s')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% objective function - single objective #1
function f1_s = myfun1(X)

% variables
d   = X(1);     % inner diameter of the hull                                (m)
t   = X(2);     % thickness of the hull                                     (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% parameters
[~, rho, mu, ~, ~, ~, ~, v, ~, ~, theta, alpha, ~, l, w] = set_Params();

% good and bad values
f1_g    = 23.6389;
f1_b    = 25.1932;

% Force of Drag (N)
f1	= eval_f1(rho, mu, v, theta, alpha, l, w, d, t, L);  

% scalled values
f1_s = (f1 - f1_g)/(f1_b - f1_g);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% objective function - single objective #2
function f2_s = myfun2(X)

% variables
d   = X(1);     % inner diameter of the hull                                (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% params
[~, ~, ~, ~, ~, ~, ~, ~, ~, ~, theta, ~, ~, ~, ~] = set_Params();

% good and bad values
f2_b    = 0.0353;
f2_g    = 0.1242;

% Internal Volume (m^3)
f2	= eval_f2(theta, d, L); 

% scalled values
f2_s = (f2 - f2_g)/(f2_b - f2_g);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% objective function - MULTIOBJECTIVE
function Lq = myfunLq(X, q, f1_0, f2_0)

% variables
d   = X(1);     % inner diameter of the hull                                (m)
t   = X(2);     % thickness of the hull                                     (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% parameters
[~, rho, mu_w, ~, ~, ~, ~, ...
    v, ~, ~, theta, alpha, ~, l, w] = set_Params();

% Force of Drag (N)
f1	= eval_f1(rho, mu_w, v, theta, alpha, l, w, d, t, L);  

% Internal Volume (m^3)
f2	= eval_f2(theta, d, L);

% good and bad values
f1_g    = 23.6389;
f1_b    = 25.1932;
f2_b    = 0.0353;
f2_g    = 0.1242;

% scalled values
f1_s = (f1 - f1_g)/(f1_b - f1_g);
f2_s = (f2 - f2_g)/(f2_b - f2_g);

% L-q norm
if isinf(q)
    Lq = max(f1_s - f1_0, f2_s - f2_0);
else
    Lq = (abs(f1_s - f1_0)^q + abs(f2_s - f2_0)^q)^(1/q);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% nonlinear constraint function
function [c, ceq] = mycon(X) 

% variables
d   = X(1);     % inner diameter of the hull                                (m)
t   = X(2);     % thickness of the hull                                     (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% parameters
[g, rho, ~, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, ~, theta, alpha, tfins, l, w] = set_Params();

% limits
[~, ~, ~, ~, ~, ~, W_lim, FS] = set_Lims();

% g7  = upper bound weight constraint
% g8 = lower bound bouyancy constraint
% g9 = upper bound bouyancy constraint
% g10 = upper bound stress constraint
[g7, g8, g9, g10] = eval_g710(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, theta, alpha, tfins, l, w,...
    d, t, L,  ...
    W_lim, FS);

% inequality constraints
c   = [g7, g8, g9, g10];

% equality constraints
ceq = [];

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END
