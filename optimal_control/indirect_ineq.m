clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Indirect method for inequality constriant
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global rho rho0 tol t0 tf x0 xf M1 M2
rho = 8;
rho0 = 10;
tol = 1e-4;
omega = 4;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

% Initial and final states
t0 = 0;
tf = 0.25;
x0 = [-rho0; 0; 0; 0; 0; pi];
xf = [0; -rho0; 0; 0; 0; pi];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------------ Optimization ----------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% X()
function J = obj_func(X)
global rho rho0 tol t0 tf x0 xf M1 M2
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------- ODE - dynamics equations -----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%
% Unconstrained arc
%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = odefun_unc(t, y)
global rho
dydt = zeros(12, 1);
% Constant
omega = 4;                                  % angular velocity, 4 rad/h

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

% Lagrange multiplier
r = y(1:3);
v = y(4:6);
lambda13 = y(7:9);
lambda46 = y(10:12);
mu = 0;

dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v + 2 * mu * r - lambda46;
dydt(7:9) = (4 * mu) * M1 * r - M1 * lambda46 - (2 * mu) * M2 * v + (4 * mu^2) * r - (2 * mu) * lambda46;
dydt(10:12) = M2 * lambda46 - lambda13;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Constrained arc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = odefun_con(t, y)
global rho
dydt = zeros(12, 1);
% Constant
omega = 4;                                  % angular velocity, 4 rad/h

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

% Lagrange multiplier
r = y(1:3);
v = y(4:6);
lambda13 = y(7:9);
lambda46 = y(10:12);
mu = 1 / (2 * (r' * r)) * (r' * lambda46 - v' * v - r' * M1 * r - r' * M2 * v);

dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v + 2 * mu * r - lambda46;
dydt(7:9) = (4 * mu) * M1 * r - M1 * lambda46 - (2 * mu) * M2 * v + (4 * mu^2) * r - (2 * mu) * lambda46;
dydt(10:12) = M2 * lambda46 - lambda13 + (4 * mu) * v - (2 * mu) * M2 * r;
end