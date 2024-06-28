%-------------------------------------------------------------------%
% Spherically Inequality Constraint                                 %
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
path_gpops = "data/gpops_ineq_data.mat";
load(path_gpops);
J_gpops = J;
x_gpops = x;
y_gpops = y;
z_gpops = z;
[theta_gpops, phi_gpops, R_gpops] = cart2sph(x_gpops, y_gpops, z_gpops);
r_gpops = r;
u1_gpops = u1;
u2_gpops = u2;
u3_gpops = u3;
u_gpops = u;
t_gpops = t;
tSolve_gpops = tSolve;
lambda_gpops = costate;
mu_gpops = mu;
H_gpops = H;

% Indirect - Lagrange Multiplier Data
path_lag = "data/indirect_ineq_data.mat";
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
H_lag = H;
eta1_lag = eta1;
eta2_lag = eta2;

%-------------------------------------------------------------------%
%-------------------------- Error Analysis -------------------------%
%-------------------------------------------------------------------%
rho_ineq = 8.8;
res_gpops = min(r_gpops) - rho_ineq;
res_lag = min(r_lag) - rho_ineq;


%-------------------------------------------------------------------%
%------------------------------- Plot ------------------------------%
%-------------------------------------------------------------------%
num_str = '_8_8';
%------------------------------- State -----------------------------%
f=figure;
plot(t_lag, r_lag, 'LineWidth', 1.5);hold on
plot(t_gpops, R_gpops, 'LineWidth', 1.5, 'LineStyle', '--');
legend('Pontryagin', 'GPOPS-II');
title('State - pos');
saveas(f, ['fig/state_ineq' num_str],'fig');


%------------------------------ Costate ----------------------------%
f=figure;

subplot(2, 3, 1)
plot(t_lag, lambda_lag(1, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_gpops, lambda_gpops(1, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'GPOPS-II')
title('\lambda_1')

subplot(2, 3, 2)
plot(t_lag, lambda_lag(2, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_gpops, lambda_gpops(2, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'GPOPS-II')
title('\lambda_2')

subplot(2, 3, 3)
plot(t_lag, lambda_lag(3, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_gpops, lambda_gpops(3, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'GPOPS-II')
title('\lambda_3')

subplot(2, 3, 4)
plot(t_lag, lambda_lag(4, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_gpops, lambda_gpops(4, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'GPOPS-II')
title('\lambda_4')

subplot(2, 3, 5)
plot(t_lag, lambda_lag(5, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_gpops, lambda_gpops(5, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'GPOPS-II')
title('\lambda_5')

subplot(2, 3, 6)
plot(t_lag, lambda_lag(6, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_gpops, lambda_gpops(6, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'GPOPS-II')
title('\lambda_6')

saveas(f, ['fig/costate_ineq' num_str],'fig');


%------------------------------ Control ----------------------------%
f=figure;
plot(t_lag, u1_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u2_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u3_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_gpops, u1_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(t_gpops, u2_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(t_gpops, u3_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('control1 - Pontryagin', 'control2 - Pontryagin', ...
       'control3 - Pontryagin', ...
       'control1 - GPOPS-II', 'control2 - GPOPS-II', ...
       'control3 - GPOPS-II');
title('Control');
saveas(f, ['fig/control_ineq' num_str],'fig');


%------------------------ Norm of Control --------------------------%
f=figure;
plot(t_lag, u_lag, 'LineWidth', 1.5);hold on
plot(t_gpops, u_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'GPOPS-II')
title('Control - Norm');
saveas(f, ['fig/controlNorm_ineq' num_str],'fig');


%--------------------------- Multiplier ----------------------------%
f=figure;
plot(t_lag, mu_lag, 'LineWidth', 1.5); hold on
plot(t_gpops, mu_gpops, 'LineWidth', 1.5, 'LineStyle', '--'); hold on
axis equal
title('Multiplier \mu');
saveas(f, ['fig/mu_ineq' num_str],'fig');
    
f=figure;
plot(t_lag, eta1_lag, 'LineWidth', 1.5);hold on
title('Multiplier \eta');
saveas(f, ['fig/eta_ineq' num_str], 'fig');


%---------------------------- Hamilton -----------------------------%
f=figure;
plot(t_lag, H_lag, 'LineWidth', 1.5); hold on
plot(t_gpops, H_gpops, 'LineWidth', 1.5, 'LineStyle', '--'); hold on
title('Hamilton');
saveas(f, ['fig/H_ineq' num_str], 'fig');

%--------------------------- Trajectory ----------------------------%
f=figure;  
plot3(x_lag, y_lag, z_lag, 'LineWidth', 1.5);hold on
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

legend('Pontryagin', 'GPOPS-II')
title('Trajectory');
saveas(f, ['fig/trajectory_ineq' num_str],'fig');
