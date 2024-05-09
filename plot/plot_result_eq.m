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
plot(t_gpops, r_gpops, 'LineWidth', 1.5, 'LineStyle', '--');
legend('Indirect Method', 'GPOPS-II');
title('State - pos');
axis equal
saveas(f, 'fig/state_eq','fig');

%------------------------------ Costate ----------------------------%
costate = lambda;
f=figure;
plot(t, costate(1, :), 'LineWidth', 1.5);hold on
plot(t, costate(2, :), 'LineWidth', 1.5);hold on
plot(t, costate(3, :), 'LineWidth', 1.5);hold on
plot(t, costate(4, :), 'LineWidth', 1.5);hold on
plot(t, costate(5, :), 'LineWidth', 1.5);hold on
plot(t, costate(6, :), 'LineWidth', 1.5);
legend('costate1', 'costate2', 'costate3', ...
       'costate4', 'costate5', 'costate6');
title('Indirect Method - Costate');
saveas(f, 'fig/costate_eq','fig');


%------------------------------ Control ----------------------------%
f=figure;
plot(t_lag, u1_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u2_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u3_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u1_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_lag, u2_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_lag, u3_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, u1_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(t_gpops, u2_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(t_gpops, u3_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('control1 - Indirect (Pontryagin)', ...
       'control2 - Indirect (Pontryagin)', ...
       'control3 - Indirect (Pontryagin)', ...
       'control1 - Indirect (Penalty)', ...
       'control2 - Indirect (Penalty)', ...
       'control3 - Indirect (Penalty)', ...
       'control1 - GPOPS-II', ...
       'control2 - GPOPS-II', ...
       'control3 - GPOPS-II');
title('Control');
saveas(f, 'fig/control_eq','fig');

%------------------------ Norm of Control --------------------------%
f=figure;
plot(t_lag, u_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_pen, u_pen, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, u_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Indirect (Pontryagin)', 'Indirect (Penalty)', 'GPOPS-II');
title('Control Norm');
saveas(f, 'fig/controlNorm_eq','fig');

%--------------------------- Multiplier ----------------------------%
f=figure;
plot(t, mu_lag, 'LineWidth', 1.5);
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
