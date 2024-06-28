%-------------------------------------------------------------------%
% Spherically constraint                                            %
% Plot Module                                                       %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
clc;clear

%-------------------------------------------------------------------%
%---------------------------- Load Data ----------------------------%
%-------------------------------------------------------------------%

% CubicFitRot Data
path_cubic = "data/cubicfitrot_data.mat";
load(path_cubic);

% GPOPS-II Data
path_gpops = "data/gpops_eq_data.mat";
load(path_gpops);
J_gpops = J;
x_gpops = x;
y_gpops = y;
z_gpops = z;
r_gpops = r;
u1_gpops = u1;
u2_gpops = u2;
u3_gpops = u3;
u_gpops = u;
t_gpops = t;
tSolve_gpops = tSolve;

% Indirect - Pontryagin Data
path_lag = "data/indirect_eq_data.mat";
load(path_lag);
J_lag = J;
x_lag = x;
y_lag = y;
z_lag = z;
r_lag = r;
u1_lag = u1;
u2_lag = u2;
u3_lag = u3;
u_lag = u;
t_lag = t;
tSolve_lag = tSolve;
lambda_lag = lambda;
mu_lag = mu;

% Indirect - Penalty Data
path_pen = "data/indirect_eq_penalty_data.mat";
load(path_pen);
J_pen = J;
x_pen = x;
y_pen = y;
z_pen = z;
r_pen = r;
u1_pen = u1;
u2_pen = u2;
u3_pen = u3;
u_pen = u;
t_pen = t;
tSolve_pen = tSolve;
lambda_pen = lambda;
mu_pen = mu;


%-------------------------------------------------------------------%
%-------------------------- Error Analysis -------------------------%
%-------------------------------------------------------------------%
rho_ineq = 10;
res_gpops = min(r_gpops) - rho_ineq;
res_lag = min(r_lag) - rho_ineq;


%-------------------------------------------------------------------%
%------------------------------- Plot ------------------------------%
%-------------------------------------------------------------------%

%------------------------------- State -----------------------------%
f=figure;
plot(t_lag, r_lag, 'LineWidth', 1.5);hold on
plot(t_pen, r_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, r_gpops, 'LineWidth', 1.5, 'LineStyle', '--');
axis equal
legend('Pontryagin', 'Penalty', 'GPOPS-II');
title('State');
saveas(f, 'fig/state_eq','fig');

%------------------------------ Costate ----------------------------%
f=figure;
subplot(2, 3, 1)
plot(t_lag, lambda_lag(1, :), 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, lambda_pen(1, :), 'LineWidth', 1.5, 'LineStyle', '-.');hold on
legend('Pontryagin', 'Penalty')
title('\lambda_1')

subplot(2, 3, 2)
plot(t_lag, lambda_lag(2, :), 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, lambda_pen(2, :), 'LineWidth', 1.5, 'LineStyle', '-.');hold on
legend('Pontryagin', 'Penalty')
title('\lambda_2')

subplot(2, 3, 3)
plot(t_lag, lambda_lag(3, :), 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, lambda_pen(3, :), 'LineWidth', 1.5, 'LineStyle', '-.');hold on
legend('Pontryagin', 'Penalty')
title('\lambda_3')
saveas(f, 'fig/costate_eq','fig');

subplot(2, 3, 4)
plot(t_lag, lambda_lag(4, :), 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, lambda_pen(4, :), 'LineWidth', 1.5, 'LineStyle', '-.');hold on
legend('Pontryagin', 'Penalty')
title('\lambda_4')

subplot(2, 3, 5)
plot(t_lag, lambda_lag(5, :), 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, lambda_pen(5, :), 'LineWidth', 1.5, 'LineStyle', '-.');hold on
legend('Pontryagin', 'Penalty')
title('\lambda_5')

subplot(2, 3, 6)
plot(t_lag, lambda_lag(6, :), 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, lambda_pen(6, :), 'LineWidth', 1.5, 'LineStyle', '-.');hold on
legend('Pontryagin', 'Penalty')
title('\lambda_6')
saveas(f, 'fig/costate_eq','fig');

%------------------------------ Control ----------------------------%
f1=figure;
plot(t_lag, u1_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, u1_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, u1_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(tcubic, u1cubic, 'LineWidth', 1.5, 'LineStyle', ':');hold on
legend('Pontryagin', 'Penalty', 'GPOPS-II', 'CubicFitRot')
title('Control-x');

f2=figure;
plot(t_lag, u2_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, u2_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, u2_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(tcubic, u2cubic, 'LineWidth', 1.5, 'LineStyle', ':');hold on
legend('Pontryagin', 'Penalty', 'GPOPS-II', 'CubicFitRot')
title('Control-y');

f3=figure;
plot(t_lag, u3_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, u3_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, u3_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(tcubic, u3cubic, 'LineWidth', 1.5, 'LineStyle', ':');hold on
legend('Pontryagin', 'Penalty', 'GPOPS-II', 'CubicFitRot')
title('Control-z');


saveas(f1, 'fig/control1_eq','fig');
saveas(f2, 'fig/control2_eq','fig');
saveas(f3, 'fig/control3_eq','fig');

%------------------------ Norm of Control --------------------------%
f=figure;
plot(t_lag, u_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, u_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, u_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(tcubic, ucubic, 'LineWidth', 1.5, 'LineStyle', ':');hold on
legend('Indirect (Pontryagin)', 'Indirect (Penalty)', 'GPOPS-II', 'CubicFitRot');
title('Control Norm');
saveas(f, 'fig/controlNorm_eq','fig');

%--------------------------- Multiplier ----------------------------%
f=figure;
plot(t_lag, mu_lag, 'LineWidth', 1.5);
title('mu');
saveas(f, 'fig/mu_eq','fig');

%--------------------------- Trajectory ----------------------------%
f=figure;  
plot3(x_lag, y_lag, z_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot3(x_pen, y_pen, z_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot3(x_gpops, y_gpops, z_gpops, ...
      'LineWidth', 1.5, 'LineStyle', '--');hold on

plot3(0, 0, 0, 'k*', 'LineWidth', 3);hold on
text(0, 0, 0, 'Chief');hold on
plot3(x_lag(1), y_lag(1), z_lag(1), 'g*', 'LineWidth', 1.5);hold on
text(x_lag(1), y_lag(1), z_lag(1), 'Departure');hold on
plot3(x_lag(end), y_lag(end), z_lag(end), 'c*', 'LineWidth', 1.5);hold on
text(x_lag(end), y_lag(end), z_lag(end), 'Arrival');hold on

rb = rho;
index = 1:1000:size(u, 2);
[X, Y, Z] = sphere;
X2 = X * rb;
Y2 = Y * rb;
Z2 = Z * rb;
surf(X2, Y2, Z2,  'FaceAlpha', 0.2, 'EdgeColor', 'texturemap'); hold on
colormap(gca, 'gray')
axis equal

legend('Indirect (Pontryagin)', 'Indirect (Penalty)', 'GPOPS-II');
title('Trajectory');
saveas(f, 'fig/trajectory_eq','fig');
