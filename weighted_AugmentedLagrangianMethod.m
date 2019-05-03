
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "weighted_AugmentedLagrangianMethod"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%       Uses a custom made Augmented Lagrangian Method 
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
%       And the results are saved in 'ALM_results.mat'
%       Currently, saving of the plots are commented out
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% set up

close all
clear
clc

% limits
[g, rho, mu_w, rho_load, rho_fins, rho_hull, Sy_hull, v, depth, T, theta, alpha, tfins, l, w] = set_Params();
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

% good and bad values
f1_g = 23.6389;
f1_b = 25.1932;
f2_b = 0.0353;
f2_g = 0.1242;

% variable bounds
LB = [d_L, t_L, L_L]; % Variable Lower Bounds
UB = [d_U, t_U, L_U]; % Variable Upper Bounds

% initial guess
d0 = d_L + (d_U - d_L)*(0);
t0 = t_L + (t_U - t_L)*(0);
L0 = L_L + (L_U - L_L)*(0);
X0 = [d0, t0, L0];

%% algorithm paramaters

% initial mu
mu = zeros(10,1);

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
        [X_min, f_min] = fminunc(@(X)fun(X, rp, mu, w1, w2), X0);
        % constraint check& parameter update
         [INFESIBILITY, rp, mu] = conup(X_min, rp, mu, gamma);

        % update penalty
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
        f1(j,:) = eval_f1(rho, mu_w, v, theta, alpha, l, w, d, t, L);  

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
% saveas(gcf,'ALM_NormalizedCriterionSpace.jpg')

% normalized criterion vs weight
figure; hold on
plot(W1, f1_s, '.')
plot(W1, f2_s, '.')
xlabel('W1')
legend('f1 - drag', 'f2 - volume')
title('Normalized Criterion Space vs Weight')
axis([0,1,0,1])
% saveas(gcf,'ALM_NormalizedCriterionSpace_WeightedTrend.jpg')

%% save results
save('ALM_results.mat', 'X', 'f1', 'f2', 'f1_s', 'f2_s', 'W1', 'W2', 'Lq1', 'Lq2', 'Lqinf')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% objective function
function A = fun(X, rp, mu, w1, w2)

% variables
d   = X(1);     % inner diameter of the hull                                (m)
t   = X(2);     % thickness of the hull                                     (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% split mu
mu1 = mu(1);
mu2 = mu(2);
mu3 = mu(3);
mu4 = mu(4);
mu5 = mu(5);
mu6 = mu(6);
mu7 = mu(7);
mu8 = mu(8);
mu9 = mu(9);
mu10 = mu(10);

% parameters
[g, rho, mu_w, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, ~, theta, alpha, tfins, l, w] = set_Params();

[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

% Force of Drag (N)
f1	= eval_f1(rho, mu_w, v, theta, alpha, l, w, d, t, L);  

% Internal Volume (m^3)
f2	= eval_f2(theta, d, L);

% good and bad values
f1_g = 23.6389;
f1_b = 25.1932;
f2_b = 0.0353;
f2_g = 0.1242;

% scalled values
f1_s = (f1 - f1_g)/(f1_b - f1_g);
f2_s = (f2 - f2_g)/(f2_b - f2_g);

% weighted sum
f = w1*f1_s + w2*f2_s;

% calc constraints
[g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = eval_gALL(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, theta, alpha, tfins, l, w, ...
    d, t, L, ...
    d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS);

% calc psi
psi1 = max(-mu1/(2*rp), g1);
psi2 = max(-mu2/(2*rp), g2);
psi3 = max(-mu3/(2*rp), g3);
psi4 = max(-mu4/(2*rp), g4);
psi5 = max(-mu5/(2*rp), g5);
psi6 = max(-mu6/(2*rp), g6);
psi7 = max(-mu7/(2*rp), g7);
psi8 = max(-mu8/(2*rp), g8);
psi9 = max(-mu9/(2*rp), g9);
psi10 = max(-mu10/(2*rp), g10);

% calc infesibility
INFESIBILITY = sumsqr([psi1,psi2,psi3,psi4,psi5,psi6,psi7,psi8,psi9,psi10]);

% calc lagrangian
L = f + ...
    mu1*psi1 + mu2*psi2 + mu3*psi3 + mu4*psi4 + mu5*psi5 + ...
    mu6*psi6 + mu7*psi7 + mu8*psi8 + mu9*psi9 + mu10*psi10;

% calc augmented lagrangian
A = L + rp*INFESIBILITY;

end

%% constraints reformulated as Infesibility
function [INFESIBILITY, rp, mu] = conup(X, rp, mu, gamma)

% variables
d   = X(1);     % inner diameter of the hull                                (m)
t   = X(2);     % thickness of the hull                                     (m)
L   = X(3);     % length of the cylindrical section of the hull             (m)

% split mu
mu1 = mu(1);
mu2 = mu(2);
mu3 = mu(3);
mu4 = mu(4);
mu5 = mu(5);
mu6 = mu(6);
mu7 = mu(7);
mu8 = mu(8);
mu9 = mu(9);
mu10 = mu(10);

% parameters
[g, rho, ~, ...
    rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, ~, theta, alpha, tfins, l, w] = set_Params();

% limits
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

% calc constraints
[g1, g2, g3, g4, g5, g6, g7, g8, g9, g10] = eval_gALL(...
    g, rho, rho_load, rho_fins, rho_hull, Sy_hull, ...
    v, depth, theta, alpha, tfins, l, w, ...
    d, t, L, ...
    d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS);

% calc psi
psi1 = max(-mu1/(2*rp), g1);
psi2 = max(-mu2/(2*rp), g2);
psi3 = max(-mu3/(2*rp), g3);
psi4 = max(-mu4/(2*rp), g4);
psi5 = max(-mu5/(2*rp), g5);
psi6 = max(-mu6/(2*rp), g6);
psi7 = max(-mu7/(2*rp), g7);
psi8 = max(-mu8/(2*rp), g8);
psi9 = max(-mu9/(2*rp), g9);
psi10 = max(-mu10/(2*rp), g10);

% calc infesibility
INFESIBILITY = sumsqr([psi1,psi2,psi3,psi4,psi5,psi6,psi7,psi8,psi9,psi10]);

% update mu's
mu1 = mu1 + 2*rp*psi1;
mu2 = mu2 + 2*rp*psi2;
mu3 = mu3 + 2*rp*psi3;
mu4 = mu4 + 2*rp*psi4;
mu5 = mu5 + 2*rp*psi5;
mu6 = mu6 + 2*rp*psi6;
mu7 = mu7 + 2*rp*psi7;
mu8 = mu8 + 2*rp*psi8;
mu9 = mu9 + 2*rp*psi9;
mu10 = mu10 + 2*rp*psi10;

% reform mu
mu = [mu1, mu2, mu3, mu4, mu5, mu6, mu7, mu8, mu9, mu10];

% update penalty
rp = gamma*rp;

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% END
