%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% script PlotGridEval
% ENME 610 - Engineering Optimization
% University of Maryland, College Park
% Group 1: David Smart, Luke Travisiano, Jason Morin
% AUV Optimization
%
%% Description:
%       Generates plots of the feasible domain in the design space and 
%       normalized criterion space where color represents how good the
%       point is according to a certain metric: f1, f2, L1, L2, Linf
%
%% Instructions:
%       (check that 'GRID_results.mat' has been created and is in the 
%       current folder)
%       Just hit "Run".
%       It will generate the plots...but it will take a while...
%       Currently, saveing of the plots and rotating gifs is commented out.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% clean up
close all
clear
clc

[d_L, d_U, t_L, t_U, L_L, L_U, ~, ~] = set_Lims();

%% GRID
GRID_res = load('GRID_results.mat');
X = GRID_res.X;
f1 = GRID_res.f1;
f2 = GRID_res.f2;
f1_s = GRID_res.f1_s;
f2_s = GRID_res.f2_s;
Lq1 = GRID_res.Lq1;
Lq2 = GRID_res.Lq2;
LqInf = GRID_res.Lqinf;

f1_g    = min(f1);
f1_b    = max(f1);
f2_b    = min(f2);
f2_g    = max(f2);

%% design space -----------------------------------------------------------
figure;
plot3(X(:,1), X(:,2), X(:,3), 'color', [0.5, 0.5, 0.5])
axis([d_L, d_U, t_L, t_U, L_L, L_U])
xlabel('d')
ylabel('t')
zlabel('L')
title('Design Space')
% saveas(gcf, 'GridEval_DesignSpace.jpg')

% % make gif
% make_animated_gif('clear')
% for i = 1:360
%     view(-45 + i,15)
%     make_animated_gif('snap', gcf)
% end
% make_animated_gif('write','GridEval_DesignSpace', 0.1)

%% in terms of f1_s
figure; hold on
for i = 1:max(size(f1_s))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*f1_s(i); % calculate color interpolated between green (0) and red (1)
    plot3(X(i,1), X(i,2),X(i,3), '.','color',c)
end
xlabel('d')
ylabel('t')
zlabel('L')
title('Design Space - f1')
axis([d_L, d_U, t_L, t_U, L_L, L_U])
view(-45,15)

% mark the best value
[~, idx] = min(f1_s);
X_f1 = X(idx, :);
plot3(X_f1(1), X_f1(2),X_f1(3),'*m')
% saveas(gcf,'GridEval_DesignSpace_f1.jpg')

% % make gif
% make_animated_gif('clear')
% for i = 1:360
%     view(-45 + i,15)
%     make_animated_gif('snap', gcf)
% end
% make_animated_gif('write','GridEval_DesignSpace_f1', 0.1)

%% in terms of f2_s
figure; hold on
for i = 1:max(size(f2_s))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*f2_s(i); % calculate color interpolated between green (0) and red (1)
    plot3(X(i,1), X(i,2),X(i,3), '.','color',c)
end
xlabel('d')
ylabel('t')
zlabel('L')
title('Design Space - f2')
axis([d_L, d_U, t_L, t_U, L_L, L_U])
view(-45,15)

% mark the best value
[~, idx] = min(f2_s);
X_f2 = X(idx, :);
plot3(X_f2(1), X_f2(2),X_f2(3),'*m')
% saveas(gcf,'GridEval_DesignSpace_f2.jpg')

% % make gif
% make_animated_gif('clear')
% for i = 1:360
%     view(-45 + i,15)
%     make_animated_gif('snap', gcf)
% end
% make_animated_gif('write','GridEval_DesignSpace_f2', 0.1)

%% in terms of Lq-1
Lq1 = (Lq1 - min(Lq1))/(max(Lq1) - min(Lq1)); % normalize to be between 0 and 1
figure; hold on
for i = 1:max(size(Lq1))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*Lq1(i); % calculate color interpolated between green (0) and red (1)
    plot3(X(i,1), X(i,2),X(i,3), '.','color',c)
end
xlabel('d')
ylabel('t')
zlabel('L')
title('Design Space - Lq1')
axis([d_L, d_U, t_L, t_U, L_L, L_U])
view(-45,15)

% mark the best value
[~, idx] = min(Lq1);
X_Lq1 = X(idx, :);
plot3(X_Lq1(1), X_Lq1(2),X_Lq1(3),'*m')
% saveas(gcf,'GridEval_DesignSpace_Lq1.jpg')

% % make gif
% make_animated_gif('clear')
% for i = 1:360
%     view(-45 + i,15)
%     make_animated_gif('snap', gcf)
% end
% make_animated_gif('write','GridEval_DesignSpace_Lq1', 0.1)

%% in terms of Lq-2
Lq2 = (Lq2 - min(Lq2))/(max(Lq2) - min(Lq2)); % normalize to be between 0 and 1
figure; hold on
for i = 1:max(size(Lq2))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*Lq2(i); % calculate color interpolated between green (0) and red (1)
    plot3(X(i,1), X(i,2),X(i,3), '.','color',c)
end
xlabel('d')
ylabel('t')
zlabel('L')
title('Design Space - Lq2')
axis([d_L, d_U, t_L, t_U, L_L, L_U])
view(-45,15)

% mark the best value
[~, idx] = min(Lq2);
X_Lq2 = X(idx, :);
plot3(X_Lq2(1), X_Lq2(2),X_Lq2(3),'*m')
% saveas(gcf,'GridEval_DesignSpace_Lq2.jpg')

% % make gif
% make_animated_gif('clear')
% for i = 1:360
%     view(-45 + i,15)
%     make_animated_gif('snap', gcf)
% end
% make_animated_gif('write','GridEval_DesignSpace_Lq21', 0.1)

%% in terms of Lq-inf
LqInf = (LqInf - min(LqInf))/(max(LqInf) - min(LqInf)); % normalize to be between 0 and 1
figure; hold on
for i = 1:max(size(LqInf))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*LqInf(i); % calculate color interpolated between green (0) and red (1)
    plot3(X(i,1), X(i,2),X(i,3), '.','color',c)
end
xlabel('d')
ylabel('t')
zlabel('L')
title('Design Space - LqInf')
axis([d_L, d_U, t_L, t_U, L_L, L_U])
view(-45,15)

% mark the best value
[~, idx] = min(LqInf);
X_Lqinf = X(idx, :);
plot3(X_Lqinf(1), X_Lqinf(2),X_Lqinf(3),'*m')
% saveas(gcf,'GridEval_DesignSpace_LqInf.jpg')

% % make gif
% make_animated_gif('clear')
% for i = 1:360
%     view(-45 + i,15)
%     make_animated_gif('snap', gcf)
% end
% make_animated_gif('write','GridEval_DesignSpace_LqInf', 0.1)

%% normalized corterion space ---------------------------------------------
figure;
plot(f1_s, f2_s, 'color', [0.5, 0.5, 0.5])
xlabel('f1 - drag')
ylabel('f2 - volume')
title('Normalized Criterion Space')
axis([0, 1, 0, 1])
% saveas(gcf, 'GridEval_NormalizedCriterionSpace.jpg')

%% in terms of f1_s
figure; hold on
for i = 1:max(size(f1_s))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*f1_s(i); % calculate color interpolated between green (0) and red (1)
    plot(f1_s(i), f2_s(i),X(i), '.','color',c)
end
xlabel('f1 - drag')
ylabel('f2 - volume')
title('Normalized Criterion Space - f1')
axis([0, 1, 0, 1])
view(-45,15)

% mark the best value
plot(1, 0,'*m')
% saveas(gcf,'GridEval_NormalizedCriterionSpace_f1.jpg')

%% in terms of f2_s
figure; hold on
for i = 1:max(size(f2_s))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*f2_s(i); % calculate color interpolated between green (0) and red (1)
    plot(f1_s(i), f2_s(i),X(i), '.','color',c) 
end
xlabel('f1 - drag')
ylabel('f2 - volume')
title('Normalized Criterion Space - f2')
axis([0, 1, 0, 1])
view(-45,15)

% mark the best value
plot(0, 1,'*m')
% saveas(gcf,'GridEval_NormalizedCriterionSpace_f2.jpg')

%% in terms of Lq-1
Lq1 = (Lq1 - min(Lq1))/(max(Lq1) - min(Lq1)); % normalize to be between 0 and 1
figure; hold on
for i = 1:max(size(Lq1))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*Lq1(i); % calculate color interpolated between green (0) and red (1)
    plot(f1_s(i), f2_s(i),X(i), '.','color',c)
end
xlabel('f1 - drag')
ylabel('f2 - volume')
title('Normalized Criterion Space - Lq1')
axis([0, 1, 0, 1])
view(-45,15)

% mark the best value
[~, idx] = min(Lq1);
f1_Lq1 = f1(idx);
f2_Lq1 = f2(idx);
plot(f1_Lq1, f2_Lq1,'*m')
% saveas(gcf,'GridEval_NormalizedCriterionSpace_Lq1.jpg')

%% in terms of Lq-2
Lq2 = (Lq2 - min(Lq2))/(max(Lq2) - min(Lq2)); % normalize to be between 0 and 1
figure; hold on
for i = 1:max(size(Lq2))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*Lq2(i); % calculate color interpolated between green (0) and red (1)
    plot(f1_s(i), f2_s(i),X(i), '.','color',c)
end
xlabel('f1 - drag')
ylabel('f2 - volume')
title('Normalized Criterion Space - Lq2')
axis([0, 1, 0, 1])
view(-45,15)

% mark the best value
[~, idx] = min(Lq2);
f1_Lq1 = f1(idx);
f2_Lq1 = f2(idx);
plot(f1_Lq1, f2_Lq1,'*m')
% saveas(gcf,'GridEval_NormalizedCriterionSpace_Lq2.jpg')

%% in terms of Lq-inf
LqInf = (LqInf - min(LqInf))/(max(LqInf) - min(LqInf)); % normalize to be between 0 and 1
figure; hold on
for i = 1:max(size(LqInf))
    c = [0,1,0] + ([1,0,0] - [0,1,0])*LqInf(i); % calculate color interpolated between green (0) and red (1)
    plot(f1_s(i), f2_s(i),X(i), '.','color',c)
end
xlabel('f1 - drag')
ylabel('f2 - volume')
title('Normalized Criterion Space - LqInf')
axis([0, 1, 0, 1])
view(-45,15)

% mark the best value
[~, idx] = min(LqInf);
f1_Lq1 = f1(idx);
f2_Lq1 = f2(idx);
plot(f1_Lq1, f2_Lq1,'*m')
% saveas(gcf,'GridEval_NormalizedCriterionSpace_LqInf.jpg')

%% END