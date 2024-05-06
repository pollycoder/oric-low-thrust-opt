%-------------------------------------------------------------------%
% Spherically constrained - GPOPS-II                                %
% Plot Module                                                       %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
clc;clear

%-------------------------------------------------------------------%
%---------------------------- Load Data ----------------------------%
%-------------------------------------------------------------------%

% GPOPS-II Data
path_gpops = "data/gpops_data.mat";
load(path_gpops);
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

% Indirect - Lagrange Multiplier Data
path_lag = "data/indirect_ineq_data.mat";
load(path_lag);
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
plot(t_gpops, r_gpops, 'LineWidth', 1.5);
legend('Indirect Method', 'GPOPS-II');
title('State - pos');
axis equal
saveas(f, 'fig/state','fig');

%------------------------------ Costate ----------------------------%
costate = lambda;
f=figure;
plot(t, costate(1, :), 'LineWidth', 1.5);hold on
plot(t, costate(2, :), 'LineWidth', 1.5);hold on
plot(t, costate(3, :), 'LineWidth', 1.5);hold on
plot(t, costate(4, :), 'LineWidth', 1.5);hold on
plot(t, costate(5, :), 'LineWidth', 1.5);hold on
plot(t, costate(6, :), 'LineWidth', 1.5);
legend('costate1', 'costate2', 'costate3', 'costate4', 'costate5', 'costate6');
title('Indirect Method - Costate');
saveas(f, 'fig/costate','fig');


%------------------------------ Control ----------------------------%
f=figure;
plot(t_lag, u1_lag, 'LineWidth', 1.5, ...
    'LineStyle', '-', 'Color', "#A2142F");hold on
plot(t_lag, u2_lag, 'LineWidth', 1.5, ...
    'LineStyle', '-', 'Color', "#0072BD");hold on
plot(t_lag, u3_lag, 'LineWidth', 1.5, ...
    'LineStyle', '-', 'Color', "#EDB120");hold on
plot(t_gpops, u1_gpops, 'LineWidth', 1.5, ...
    'LineStyle', '--', 'Color', "#A2142F");hold on
plot(t_gpops, u2_gpops, 'LineWidth', 1.5, ...
    'LineStyle', '--', 'Color', "#0072BD");hold on
plot(t_gpops, u3_gpops, 'LineWidth', 1.5, ...
    'LineStyle', '--', 'Color', "#EDB120");hold on
legend('control1 - Indirect', 'control2 - Indirect', ...
       'control3 - Indirect', 'control1 - GPOPS-II', ...
       'control2 - GPOPS-II', 'control3 - GPOPS-II');
title('Control');
saveas(f, 'fig/control','fig');

%------------------------ Norm of Control --------------------------%
f=figure;
plot(t_lag, u_lag, 'LineWidth', 1.5, 'Color', "#0072BD");hold on
plot(t_gpops, u_gpops, 'LineWidth', 1.5, 'Color', "#A2142F");hold on
legend('control - Indirect', 'control - GPOPS-II');
title('Control - Norm');
saveas(f, 'fig/controlNorm','fig');

%--------------------------- Multiplier ----------------------------%
f=figure;
plot(t, mu_lag, 'LineWidth', 1.5);
title('mu');
saveas(f, 'fig/mu','fig');

%--------------------------- Trajectory ----------------------------%
f=figure;  
plot3(x_lag, y_lag, z_lag, 'LineWidth', 1.5, 'Color', "#A2142F");hold on
plot3(x_gpops, y_gpops, z_gpops, 'LineWidth', 1.5, 'Color', "#0072BD");hold on

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

legend('Indirect', 'GPOPS-II');
title('Trajectory');
saveas(f, 'fig/trajectory','fig');
