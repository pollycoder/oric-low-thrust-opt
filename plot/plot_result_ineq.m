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
eta1_lag = eta1;
eta2_lag = eta2;

% Indirect - Lagrange Multiplier Data - Switch Function
path_lag = "data/indirect_ineq_sf_data.mat";
load(path_lag);
J_lag_sf = J;
x_lag_sf = x;
y_lag_sf = y;
z_lag_sf = z;
r_lag_sf = r;
u1_lag_sf = u1;
u2_lag_sf = u2;
u3_lag_sf = u3;
u_lag_sf = u;
t_lag_sf = t;
tSolve_lag_sf = tSolve;
lambda_lag_sf = lambda;
mu_lag_sf = mu;
eta1_lag_sf = eta1;
eta2_lag_sf = eta2;



%-------------------------------------------------------------------%
%-------------------------- Error Analysis -------------------------%
%-------------------------------------------------------------------%
rho_ineq = 9;
res_gpops = min(r_gpops) - rho_ineq;
res_lag = min(r_lag) - rho_ineq;


%-------------------------------------------------------------------%
%------------------------------- Plot ------------------------------%
%-------------------------------------------------------------------%

%------------------------------- State -----------------------------%
f=figure;
plot(t_lag, r_lag, 'LineWidth', 1.5);hold on
plot(t_lag_sf, r_lag_sf, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, R_gpops, 'LineWidth', 1.5, 'LineStyle', '--');
legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II');
title('State - pos');
saveas(f, 'fig/state_ineq','fig');

%------------------------------ Costate ----------------------------%
f=figure;

subplot(2, 3, 1)
plot(t_lag, lambda_lag(1, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag_sf, lambda_lag_sf(1, :), ...
     'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, lambda_gpops(1, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II')
title('\lambda_1')

subplot(2, 3, 2)
plot(t_lag, lambda_lag(2, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag_sf, lambda_lag_sf(2, :), ...
     'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, lambda_gpops(2, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II')
title('\lambda_2')

subplot(2, 3, 3)
plot(t_lag, lambda_lag(3, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag_sf, lambda_lag_sf(3, :), ...
     'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, lambda_gpops(3, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II')
title('\lambda_3')

subplot(2, 3, 4)
plot(t_lag, lambda_lag(4, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag_sf, lambda_lag_sf(4, :), ...
     'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, lambda_gpops(4, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II')
title('\lambda_4')

subplot(2, 3, 5)
plot(t_lag, lambda_lag(5, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag_sf, lambda_lag_sf(5, :), ...
     'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, lambda_gpops(5, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II')
title('\lambda_5')

subplot(2, 3, 6)
plot(t_lag, lambda_lag(6, :), ...
     'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag_sf, lambda_lag_sf(6, :), ...
     'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, lambda_gpops(6, :), ...
     'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II')
title('\lambda_6')

saveas(f, 'fig/costate_ineq','fig');


%------------------------------ Control ----------------------------%
f=figure;
plot(t_lag, u1_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u2_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u3_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag_sf, u1_lag_sf, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_lag_sf, u2_lag_sf, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_lag_sf, u3_lag_sf, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, u1_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(t_gpops, u2_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
plot(t_gpops, u3_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('control1 - Pontryagin', 'control2 - Pontryagin', ...
       'control3 - Pontryagin', ...
       'control1 - Pontryagin', 'control2 - Pontryagin', ...
       'control3 - Pontryagin', ...
       'control1 - GPOPS-II', 'control2 - GPOPS-II', ...
       'control3 - GPOPS-II');
title('Control');
saveas(f, 'fig/control_ineq','fig');

%------------------------ Norm of Control --------------------------%
f=figure;
plot(t_lag, u_lag, 'LineWidth', 1.5);hold on
plot(t_lag_sf, u_lag_sf, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
plot(t_gpops, u_gpops, 'LineWidth', 1.5, 'LineStyle', '--');hold on
legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II')
title('Control - Norm');
saveas(f, 'fig/controlNorm_ineq','fig');

%--------------------------- Multiplier ----------------------------%
f=figure;
plot(t_lag, mu_lag, 'LineWidth', 1.5); hold on
plot(t_lag_sf, mu_lag_sf, 'LineWidth', 1.5);
legend('Pontryagin', 'Pontryagin - Swtich Function')
title('Multiplier \mu');
saveas(f, 'fig/mu_ineq','fig');
    
f=figure;
plot(t_lag, eta1_lag, 'LineWidth', 1.5);hold on
plot(t_lag, eta2_lag, 'LineWidth', 1.5);hold on
legend('Pontryagin - eta1', 'Pontryagin - eta2')
title('Multiplier \eta');

f=figure;
plot(t_lag_sf, eta1_lag_sf, 'LineWidth', 1.5);hold on
plot(t_lag_sf, eta2_lag_sf, 'LineWidth', 1.5);hold on
legend('Pontryagin - Swtich Function - eta1', ...
       'Pontryagin - Swtich Function - eta2')

title('Multiplier \eta - SF');
saveas(f, 'fig/eta_ineq','fig');
%%
%--------------------------- Trajectory ----------------------------%
f=figure;  
plot3(x_lag, y_lag, z_lag, 'LineWidth', 1.5);hold on
plot3(x_lag_sf, y_lag_sf, z_lag_sf, 'LineWidth', 1.5, 'LineStyle', '-.');hold on
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

legend('Pontryagin', 'Pontryagin - Switch Function', 'GPOPS-II')
title('Trajectory');
saveas(f, 'fig/trajectory_ineq','fig');
