%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "Plot_Opt_on_Feas"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%       Evaluates the inequality constraints for a full 3D cubic grid of 
%       points spanning the variable bounds. Clusters the evaluated points
%       according to the inequality constraints violated. Generates a plot
%       showing the various regions.
%           The solid green region is the feasible domain.
%           The transparent grey region violates the stress constraint. 
%           The transparent blue region violates the bouyancy constraints. 
%           The transparent red region violates the weight constraint. 
%           Note that these infeasible regions overlap.
%       Also plots the optima saved in 'fmincon_results.mat' on top.
%  
%% Instructions:
%       (Check that 'fmincon_results.mat' has been created and is in the
%       current folder...
%       ...Note: fmincon could be easily replaced to plot other optima)
%       Just hit "Run".
%           It will generate a 3D plot of the deign space showing the feasible 
%           domain and infeasible regions.
%       Hit any button to continue when ready.
%           It will add the optima points to the plot.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% clean up
close all
clear
clc

%% Variable Definitions:
%   d        inner diameter of the hull (m)
%   t        thickness of the hull (m)
%   L        length of the cylindrical section of the hull (m)

%% get all parameters:
[g, rho, mu, ...
    rho_load, rho_fins, rho_hull,  Sy_hull, ...
    v, depth, T, theta, alpha, tfins, l, w] = set_Params();

%% Variable Bounds:
[d_L, d_U, t_L, t_U, L_L, L_U, W_lim, FS] = set_Lims();

%% EVAL for all 3 variables

N = 50;

% init
W           = zeros(N,N,N);
F_B         = zeros(N,N,N);
s_max       = zeros(N,N,N);
F_D         = zeros(N,N,N);
V_i         = zeros(N,N,N);
W_load      = zeros(N,N,N);
F_L         = zeros(N,N,N);
W_T         = zeros(N,N,N);

d = linspace(d_L, d_U, N);
t = linspace(t_L, t_U, N);
L = linspace(L_L, L_U, N);

XR = []; 
XB = [];
XK = [];
XG = [];

figure(1);
hold on

for i = 1:N
    for j = 1:N
        for k = 1:N
        
            % Calculate Everything
            [W(i,j,k), F_B(i,j,k), s_max(i,j,k), F_D(i,j,k), V_i(i,j,k), W_load(i,j,k), F_L(i,j,k)] = ...
                calc_Everything(d(i), t(j), L(k));

            W_T(i,j,k) = W(i,j,k) + W_load(i,j,k);

            % CHECK CONSTRAINTS
            
            % stress - black
            if s_max(i,j,k) - (1/FS)*Sy_hull > 0 
                if isempty(XK) 
                    XK = d(i);
                    YK = t(j);
                    ZK = L(k);
                else  
                    XK(end+1,1) = d(i);
                    YK(end+1,1) = t(j);
                    ZK(end+1,1) = L(k);
                end
            end
            % bouyancy - blue
            if F_B(i,j,k) - W_T(i,j,k) - F_L(i,j,k) > 0 || W_T(i,j,k) - F_B(i,j,k) - F_L(i,j,k) > 0 
                if isempty(XB) 
                    XB = d(i);
                    YB = t(j);
                    ZB = L(k);
                else  
                    XB(end+1,1) = d(i);
                    YB(end+1,1) = t(j);
                    ZB(end+1,1) = L(k);
                end
            end
                       
            % weight - red
            if W_T(i,j,k) - W_lim > 0                
                if isempty(XR) 
                    XR = d(i);
                    YR = t(j);
                    ZR = L(k);
                else  
                    XR(end+1,1) = d(i);
                    YR(end+1,1) = t(j);
                    ZR(end+1,1) = L(k);
                end
            end
            
            % everything passed - green
            if s_max(i,j,k) - (1/FS)*Sy_hull > 0 
            elseif F_B(i,j,k) - W_T(i,j,k) - F_L(i,j,k) > 0 || W_T(i,j,k) - F_B(i,j,k) - F_L(i,j,k) > 0
            elseif W_T(i,j,k) - W_lim > 0
            else
                if isempty(XG)
                    XG = d(i);
                    YG = t(j);
                    ZG = L(k);
                else  
                    XG(end+1,1) = d(i);
                    YG(end+1,1) = t(j);
                    ZG(end+1,1) = L(k);
                end
            end
        end
    end
end

%% create shapes

% bouyancy - blue
if ~isempty(XB)
    SB = alphaShape(XB, YB, ZB);
    SB.Alpha = 0.05;
    plot(SB,'FaceColor','b','FaceAlpha',0.25,'EdgeColor','none');
end

% weight - red
if ~isempty(XR)
    SR = alphaShape(XR, YR, ZR);
    SR.Alpha = 0.05;
    plot(SR,'FaceColor','r','FaceAlpha',0.25,'EdgeColor','none');
end

% stress - black
if ~isempty(XK)
    SK = alphaShape(XK, YK, ZK);
    SK.Alpha = 0.05;
    plot(SK,'FaceColor','k','FaceAlpha',0.25,'EdgeColor','none');
end

% everything passed - green
if ~isempty(XG)
    SG = alphaShape(XG, YG, ZG);
    % SG.Alpha = 0.05;
    plot(SG,'FaceColor','g');
end
xlabel('d (m)')
ylabel('t (m)')
zlabel('L (m)')
title('Feasible Domain - Design Space')
axis([d_L d_U t_L t_U L_L L_U])
axis square
view(-45,15)
% saveas(gcf, 'FeasibleDomain_DesignSpace.jpg')

pause();

%% load in optima
ALL_res     = load('fmincon_results.mat'); % this can be changed
X           = ALL_res.X;

%% plot optima
plot3(X(:,1), X(:,2), X(:,3), 'ok', 'MarkerFaceColor', 'k')
title('Optima - Design Space')
% saveas(gcf, 'Optima_on_FeasibleDomain_DesignSpace.jpg')

% %% make gif
% make_animated_gif('clear')
% 
% for i = 1:360
%     view(-45 + i,15)
%     make_animated_gif('snap', gcf)
% end
% make_animated_gif('write','Optima_on_FeasibleDomain_DesignSpace', 0.1)

%% END