
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script "isolate_Pareto"
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%       Isolates the Pareto Frontier from a set of points
%% Instructions:
%       (Check that 'GRID_results.mat' and 'ALL_results.mat' have been
%       created and were saved in the current folder.)
%       Just hit "Run".
%       It will plot all the points in the feasible domain in grey.
%       It will plot the Pareto Frontier in green.
%       It will save the plot as 'Pareto.jpg'
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clean up
close all
clear
clc

[d_L, d_U, t_L, t_U, L_L, L_U, ~, ~] = set_Lims();

%% Grid Eval
GRID_res     = load('GRID_results.mat');
f1_s        = GRID_res.f1_s;
f2_s        = GRID_res.f2_s;

%% plot grid

% normalized corterion space
figure(3); hold on
plot(f1_s, f2_s,'.', 'color', [0.5, 0.5, 0.5])
xlabel('f1')
ylabel('f2')
title('normalized criterion space')

%% All Optima
ALL_res     = load('ALL_results.mat');
f1_s        = ALL_res.f1_s;
f2_s        = ALL_res.f2_s;

%% find the pareto points

k = 0;
for i = max(size(f1_s)):-1:1
    np = 0;
    for j = 1:i-1 % check all other values
        if f1_s(i) > f1_s(j) && f2_s(i) > f1_s(j) % if dominated by another point in terms of both objective functions
            % not pareto
            np = 1;
        end
    end
    if ~np
        % pareto
        k = k+1;
        f1_s_p(k,:) = f1_s(i,:);
        f2_s_p(k,:) = f2_s(i,:);
    end   
end

% update
f1_s = f1_s_p;
f2_s = f2_s_p;
f1_s_p = [];
f2_s_p = [];

k = 0;
for i = max(size(f1_s)):-1:1
    np = 0;
    for j = 1:i-1 % check all other values
        if f1_s(i) == f1_s(j) && f2_s(i) > f1_s(j) || f2_s(i) == f2_s(j) && f1_s(i) > f1_s(j) % if dominated by another point in terms of just one of the objectives
            % not pareto
            np = 1;
        end
    end
    if ~np
        % pareto
        k = k+1;
        f1_s_p(k,:) = f1_s(i,:);
        f2_s_p(k,:) = f2_s(i,:);
    end   
end

%% plot Pareto points
plot(f1_s_p, f2_s_p,'.g')

%% save
saveas(gcf, 'Pareto.jpg')

%%