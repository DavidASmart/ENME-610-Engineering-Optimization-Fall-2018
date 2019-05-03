
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "weighted_fimincon_script"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%       Uses the built-in MATLAB function "fmincon" in conjuntion with 
%       the Weighted Method to acheive multi-objective optimization.
%       The weight for the first objective is swept from 0 to 1.
%
%% Instructions:
%       Just hit "Run"
%       A plot of the optima acheived in the criterion space is generated.
%       A plot of the objective functions vs the weight
%       assigned to the first objective are generated.
%       And the results are saved in 'fmincon_results.mat'
%       Currently, saving of the plots are commented out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% set up

close all
clear
clc

% limits
[g, rho, mu, rho_load, rho_fins, rho_hull, Sy_hull, v, depth, T, theta, alpha, tfins, l, w] = set_Params();
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

% variable bounds
LB = [d_L, t_L, L_L]; % Variable Lower Bounds
UB = [d_U, t_U, L_U]; % Variable Upper Bounds

% initial guess
d0 = d_L + (d_U - d_L)*(1);
t0 = t_L + (t_U - t_L)*(1);
L0 = L_L + (L_U - L_L)*(1);
X0 = [d0, t0, L0];

% good and bad values
f1_g    = 23.6389;
f1_b    = 25.1932;
f2_b    = 0.0353;
f2_g    = 0.1242;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% optimization / minimization

% currently these options are pretty much default, but I put this here for
% potential uses in the future
options = optimoptions(@fmincon, ...
    'ConstraintTolerance',      1e-6, ...
    'MaxFunctionEvaluations',   2000, ...
    'MaxIterations',            1500, ...
    'StepTolerance',            1e-10);

j = 0;
for i = 0:100
    
    % multiobjective weights
    w1 = i/100;
    w2 = 1 - w1;

    % do it!
    [X_opt, ~] = fmincon(@(X)myfun(X, w1, w2), X0, [], [], [], [], LB, UB, @mycon, options);
    
    d = X_opt(1);
    t = X_opt(2);
    L = X_opt(3);
    
    % check feasibility
    [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = eval_gALL(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, theta, alpha, tfins, l, w, ...
    d, t, L, ...
    d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS);

    if g1 > 0 || g2  > 0 || g3 > 0 || g4 > 0 || g5 > 0 || g6 > 0 || g7 > 0 || g8 > 0 || g9 > 0 || g10 > 0
        
        % constraints violated
        
    else
        j = j+1;
        
        % params
        W1(j,:) = w1;
        W2(j,:) = w2;
        
        % state
        X(j,:) = X_opt;

        % Force of Drag (N)
        f1(j,:)	= eval_f1(rho, mu, v, theta, alpha, l, w, d, t, L);  

        % Internal Volume (m^3)
        f2(j,:)	= eval_f2(theta, d, L); 
    end
    
end

% update good and bad
if min(f1) < f1_g
    f1_g = min(f1);
end
if max(f1) > f1_b
    f1_b = max(f1);
end
if min(f2) < f2_b
    f2_b = min(f2);
end
if max(f2) > f2_g
    f2_g = max(f2);
end

% normalize
for i = 1:j
    % scalled values
    f1_s(i,:) = (f1(i) - f1_g)/(f1_b - f1_g);
    f2_s(i,:) = (f2(i) - f2_g)/(f2_b - f2_g);

    % multi-objective Lq-method
    Lq1(i,:) = sum([f1_s(i), f2_s(i)]);
    Lq2(i,:) = sqrt(sumsqr([f1_s(i), f2_s(i)]));
    Lqinf(i,:) = max([f1_s(i), f2_s(i)]);
end

%% plots

% normalized cirterion space
figure;
plot(f1_s, f2_s,'.b')
xlabel('f1 - drag')
ylabel('f2 - volume')
title('Normalized Criterion Space')
axis([0,1,0,1])
%saveas(gcf,'fmincon_NormalizedCriterionSpace.jpg')

% normalized criterion vs weight
figure; hold on
plot(W1, f1_s, '.')
plot(W1, f2_s, '.')
xlabel('W1')
legend('f1 - drag', 'f2 - volume')
title('normalized criterion space vs weight')
axis([0,1,0,1])
% saveas(gcf,'fmincon_NormalizedCriterionSpace_WeightedTrend.jpg')

%% save results
save('fmincon_results.mat', 'X', 'f1', 'f2', 'f1_s', 'f2_s', 'W1', 'W2', 'Lq1', 'Lq2', 'Lqinf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% objective function
function f = myfun(X, w1, w2)

% variables
d   = X(1);     % inner diameter of the hull                                (m)
t   = X(2);     % thickness of the hull                                     (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% parameters
[~, rho, mu, ~, ~, ~, ~, v, ~, ~, theta, alpha, ~, l, w] = set_Params();

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

% weighted sum
f = w1*f1_s + w2*f2_s;

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
