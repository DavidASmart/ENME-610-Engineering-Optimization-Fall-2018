
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "MultiStart_weighted_PenaltyMethod"
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
%       Multi-Start is also used to overcome the discrepency between
%       solutions depending on the start point. Many diffrent starting 
%       points are tried and the best solution is selected.
%       The weight for the first objective is swept from 0 to 1.
%       For each weight, multi-start is used.
%
%% Instructions:
%       Select the parameters for the Multi-Start generator
%       m   method of point generation
%           1 - gaussian random
%           2 - halton psudo-random
%           3 - grid
%           4 - sphere
%       N   number of points to generate
%       Nd  number of points in d (for grid only)
%       Nt  number of points in t (for grid only)
%       NL  number of points in L (for grid only)
%
%       Then, hit "Run"
%       A plot of the optima acheived in the criterion space is generated.
%       A plot of the objective functions vs the weight
%       assigned to the first objective is generated.
%       And the results are saved in 'multistart_results.mat'
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

%% algorithm paramaters

% pennalty grow rate **
gamma = 1.5;

% max pennalty **
max_rp = 1e3;

% INFESIBILITY tolerance
tol = 1e-6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% generates a set of start points
m = 3;
N = 5;
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

% % plot multi-start points
% plot3(X0(:,1), X0(:,2), X0(:,3),'.')
% axis([d_L, d_U, t_L, t_U, L_L, L_U])
% xlabel('d')
% ylabel('t')
% zlabel('L')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% iterative optimization loop

jj = 0;
for ii = 0:100
    
    % multiobjective weights
    w1 = ii/100;
    w2 = 1 - w1;

    j = 0;
    X_multi = [];
    f_multi = [];
    % multi-start
    for i = 1:size(X0,1)
        
        % check feasibility
        [g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = eval_gALL(...
        g, rho, rho_load, rho_fins, rho_hull, Sy_hull, ...
        v, depth, theta, alpha, tfins, l, w, ...
        X0(i,1), X0(i,2), X0(i,3), ...
        d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS);
        
        % penalty method must start outside the feasible region
        if g1 > 0 || g2  > 0 || g3 > 0 || g4 > 0 || g5 > 0 || g6 > 0 || g7 > 0 || g8 > 0 || g9 > 0 || g10 > 0

            % intialaize penalty method
            rp = 1;
            INFESIBILITY = 1;
            INFESIBILITY_old = 1;
            nochange = 0;

            % iterative penalty method
            while INFESIBILITY > tol && nochange == 0

                % unconstrained minimization
                [X_min, f_min] = fminunc(@(X)fun(X, rp, w1, w2), X0(i,:));
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
                X_multi(j,:) = X_min;
                f_multi(j,:) = f_min;
            end
        end
    end
    
    % if a feasible solution was found, record it
    if ~isempty(f_multi)
        % find best of multi-start for this weight
        [f_opt, idx] = min(f_multi);
        
        jj = jj + 1;
        
        % params
        W1(jj,:) = w1;
        W2(jj,:) = w2;
        
        % state
        X(jj,:) = X_multi(idx,:);
        
        d = X(jj,1);
        t = X(jj,2);
        L = X(jj,3);
        
        % Force of Drag (N)
        f1(jj,:) = eval_f1(rho, mu, v, theta, alpha, l, w, d, t, L);  

        % Internal Volume (m^3)
        f2(jj,:) = eval_f2(theta, d, L); 
    end
end

% normalize
for i = 1:jj
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
% saveas(gcf,'multistart_NormalizedCriterionSpace.jpg')

% normalized criterion vs weight
figure; hold on
plot(W1, f1_s, '.')
plot(W1, f2_s, '.')
xlabel('W1')
legend('f1 - drag', 'f2 - volume')
title('Normalized Criterion Space vs Weight')
% saveas(gcf,'multistart_NormalizedCriterionSpace_WeightedTrend.jpg')

%% save results
save('multistart_results.mat', 'X', 'f1', 'f2', 'f1_s', 'f2_s', 'W1', 'W2', 'Lq1', 'Lq2', 'Lqinf')

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
