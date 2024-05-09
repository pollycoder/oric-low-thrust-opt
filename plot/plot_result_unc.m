%-------------------------------------------------------------------%
% Unconstrained                                                     %
% Plot Module                                                       %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
clc;clear

%-------------------------------------------------------------------%
%---------------------------- Load Data ----------------------------%
%-------------------------------------------------------------------%

% Indirect - Pontryagin Data
path_lag = "data/indirect_unc_data.mat";
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


%-------------------------------------------------------------------%
%------------------------------- Plot ------------------------------%
%-------------------------------------------------------------------%

%------------------------------- State -----------------------------%
f=figure;
plot(t_lag, r_lag, 'LineWidth', 1.5);hold on
title('State');
saveas(f, 'fig/state_unc','fig');

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
title('Costate');
saveas(f, 'fig/costate_unc','fig');


%------------------------------ Control ----------------------------%
f=figure;
plot(t_lag, u1_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u2_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
plot(t_lag, u3_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
legend('control1', 'control2', 'control3');
title('Control');
saveas(f, 'fig/control_unc','fig');

%------------------------ Norm of Control --------------------------%
f=figure;
plot(t_lag, u_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
title('Control Norm');
saveas(f, 'fig/controlNorm_unc','fig');

%--------------------------- Trajectory ----------------------------%
f=figure;  
plot3(x_lag, y_lag, z_lag, 'LineWidth', 1.5, 'LineStyle', '-');hold on
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
title('Trajectory');
saveas(f, 'fig/trajectory_unc','fig');
