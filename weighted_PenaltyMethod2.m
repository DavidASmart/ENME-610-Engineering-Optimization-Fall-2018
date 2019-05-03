
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "weighted_PenaltyMethod"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%       Uses a custom made Penalty Method 
%       (utilizing fminunc for unconstrained optimization)
%       in conjuntion with the Weighted Method to acheive multi-objective 
%       optimization.
%       The weight for the first objective is swept from 0 to 1.
%
%% Instructions:
%       Set the initial position X0, and then hit "Run"
%       A plot of the optima acheived in the criterion space is generated.
%       A plot of the objective functions vs the weight
%       assigned to the first objective is generated.
%       And the results are saved in 'Penalty_results.mat'
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

% good and bad values
f1_g = 23.6389;
f1_b = 25.1932;
f2_b = 0.0353;
f2_g = 0.1242;

% variable bounds
LB = [d_L, t_L, L_L]; % Variable Lower Bounds
UB = [d_U, t_U, L_U]; % Variable Upper Bounds

% initial guess                                 *** set by user ***
d0 = d_L + (d_U - d_L)*(1/2);
t0 = t_L + (t_U - t_L)*(1/2);
L0 = L_L + (L_U - L_L)*(1/2);
X0 = [d0, t0, L0];

%% algorithm paramaters

% pennalty grow rate **
gamma = 1.5;

% max pennalty **
max_rp = 1e3;

% INFESIBILITY tolerance
tol = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterative optimization loop

j = 0;
for i = 0:100
    
    % multiobjective weights
    w1 = i/100;
    w2 = 1 - w1;
        
    % intialaize penalty method
    rp = 1;
    INFESIBILITY = 1;
    INFESIBILITY_old = 1;
    nochange = 0;

    % iterative penalty method
    while INFESIBILITY > tol && nochange == 0

        % unconstrained minimization
        [X_min, f_min] = fminunc(@(X)fun(X, rp, w1, w2), X0);
        % constraint check
        [INFESIBILITY] = con(X_min);

        % update penalty
        rp = gamma*rp;
        if rp > max_rp
            rp = max_rp;
        end

        % X0(i,:) = X_min;

        if INFESIBILITY == INFESIBILITY_old
            nochange = 1;
        end
        INFESIBILITY_old = INFESIBILITY;
    end

    % save results for all multi-starts within infeasability tolerance
    if INFESIBILITY < tol
        j = j + 1;
        
        % params
        W1(j,:) = w1;
        W2(j,:) = w2;
        
        % state
        X(j,:) = X_min;
        
        d = X(j,1);
        t = X(j,2);
        L = X(j,3);
        
        % Force of Drag (N)
        f1(j,:) = eval_f1(rho, mu, v, theta, alpha, l, w, d, t, L);  

        % Internal Volume (m^3)
        f2(j,:) = eval_f2(theta, d, L); 
    end
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
% saveas(gcf,'penalty2_NormalizedCriterionSpace.jpg')

% normalized criterion vs weight
figure; hold on
plot(W1, f1_s, '.')
plot(W1, f2_s, '.')
xlabel('W1')
legend('f1 - drag', 'f2 - volume')
title('Normalized Criterion Space vs Weight')
axis([0,1,0,1])
% saveas(gcf,'penalty2_NormalizedCriterionSpace_WeightedTrend.jpg')

%% save results
save('penalty2_results.mat', 'X', 'f1', 'f2', 'f1_s', 'f2_s', 'W1', 'W2', 'Lq1', 'Lq2', 'Lqinf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% objective function
function phi = fun(X, rp, w1, w2)

% variables
d   = X(1);     % inner diameter of the hull                                (m)
t   = X(2);     % thickness of the hull                                     (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% parameters
[g, rho, mu, ...
    rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, ~, theta, alpha, tfins, l, w] = set_Params();

% limits
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

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

% constraints
[g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = eval_gALL(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, theta, alpha, tfins, l, w, ...
    d, t, L, ...
    d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS);

INFESIBILITY = ...
    max(0, g1) + max(0, g2) + max(0, g3) + max(0, g4) + max(0, g5) + ...
    max(0, g6) + max(0, g7) + max(0, g8) + max(0, g9) + max(0, g10);

phi = f + rp*INFESIBILITY;

end

%% constraints reformulated as Infesibility
function [INFESIBILITY] = con(X)

% variables
d   = X(1);     % inner diameter of the hull                                (m)
t   = X(2);     % thickness of the hull                                     (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% parameters
[g, rho, ~, ...
    rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, ~, theta, alpha, tfins, l, w] = set_Params();

% limits
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

% constraints
[g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = eval_gALL(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, theta, alpha, tfins, l, w, ...
    d, t, L, ...
    d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS);

INFESIBILITY = sumsqr([max(0, g1),max(0, g2),max(0, g3),max(0, g4),max(0, g5),...
    max(0, g6),max(0, g7),max(0, g8),max(0, g9),max(0, g10)]);

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END
