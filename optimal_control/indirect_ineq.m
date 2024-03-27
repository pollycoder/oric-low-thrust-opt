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
t1_guess = tf / 3;
t2_guess = tf * 2 / 3;
lambda_guess1 = ones(6, 1);
lambda_guess2 = ones(6, 1);
xInit = [t1_guess; lambda_guess1; lambda_guess2; t2_guess];
[x, J] = fmincon(@obj_func, xInit, [], [], [], [], [], [], @(x)nonlcon(x));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --------------------------- Calculate and plot -------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = x(1);
t2 = x(14);

n1 = 100;
n2 = 100;
t01_span = linspace(t0, t1, n1);
t21_span = linspace(tf, t1, n2);

yt0 = [x0; x(2:7)];
ytf = [xf; x(8:13)];

[t01, y01] = ode45(@odefun_unc, t01_span, yt0);
[t21, y21] = ode45(@odefun_con, t21_span, ytf);

% State and costate
y01 = y01';
t1f = flip(t21)';
y1f = flip(y21, 1)';

% Lagrange multiplier mu (t1, tf)
r = y1f(1:3, :);
v = y1f(4:6, :);
lambda46 = y1f(10:12, :); 

mu = zeros(size(y1f, 2), 1);
for i = 1:size(y1f, 2)
    mu(i) = 1 / (2 * rho^2) * (r(:, i)' * lambda46(:, i) ...
            - v(:, i)' * v(:, i) - r(:, i)' * M1 * r(:, i) ...
            - r(:, i)' * M2 * v(:, i));
end

% Control
u01 = -y01(10:12, :);
u01_norm = vecnorm(u01);
u1f = -lambda46 + 2 * (r * mu);
u1f_norm = vecnorm(u1f);
u_norm = [u01_norm, u1f_norm];

r1f = r;
r01 = y01(1:3, :);
r1f_norm = vecnorm(r1f);
r01_norm = vecnorm(r01);
r_norm = [r01_norm, r1f_norm];
t = [t01; t1f'];

% plot
plot(t, r_norm);



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ------------------------- Optimization module --------------------------%                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% About x:                                                                
% x(1): t1                                                                
% x(2:7): lambda-guess 1 (t0)                                            
% x(8:13): lambda-guess 2 (tf) 
% x(14): t2
function J = obj_func(x)
global x0 xf t0 tf rho M2 M1
t1 = x(1);
t2 = x(14);

n1 = 100;
n2 = 100;
n3 = 100;
t01_span = linspace(t0, t1, n1);
t21_span = linspace(t2, t1, n2);
tf2_span = linspace(tf, t2, n3);

yt0 = [x0; x(2:7)];
ytf = [xf; x(8:13)];

[t01, y01] = ode45(@odefun_unc, t01_span, yt0);
[tf2, yf2] = ode45(@odefun_unc, tf2_span, ytf);

yt2 = yf2(end, :)';
[t21, y21] = ode45(@odefun_con, t21_span, yt2);


% State and costate
y01 = y01';
y12 = flip(y21, 1)';
t12 = flip(t21)';
y2f = flip(yf2, 1)';
t2f = flip(tf2)';

% Lagrange multiplier mu (t1, t2)
r = y12(1:3, :);
v = y12(4:6, :);
lambda46 = y12(10:12, :); 

mu = zeros(size(y12, 2), 1);
for i = 1:size(y12, 2)
    mu(i) = 1 / (2 * rho^2) * (r(:, i)' * lambda46(:, i) ...
            - v(:, i)' * v(:, i) - r(:, i)' * M1 * r(:, i) ...
            - r(:, i)' * M2 * v(:, i));
end

% Control
u01 = -y01(10:12, :);
u01_norm = vecnorm(u01);
u12 = -lambda46 + 2 * (r * mu);
u12_norm = vecnorm(u12);
u2f = -y2f(10:12, :);
u2f_norm = vecnorm(u2f);
u_norm = [u01_norm, u12_norm, u2f_norm];

% Integration
t_int = [t01; t12'; t2f'];
dJ = 0.5 .* u_norm .* u_norm;
J = trapz(t_int, dJ);
end


% Non-linear constraint
function [c, ceq] = nonlcon(x)
global M1 M2 rho x0 xf t0 tf
c = [];                                             % Inequality constraint

t1 = x(1);

n1 = 100;
n2 = 100;
t01_span = linspace(t0, t1, n1);
t21_span = linspace(tf, t1, n2);

yt0 = [x0; x(2:7)];
ytf = [xf; x(8:13)];

[t01, y01] = ode45(@odefun_unc, t01_span, yt0);
[t21, y21] = ode45(@odefun_con, t21_span, ytf);

rv_t1p = y21(end, 1:6)';
lambda_t1m = y01(end, 7:12)';
lambda_t1p = y21(end, 7:12)';

r_t1p = rv_t1p(1:3);
v_t1p = rv_t1p(4:6);

% Lagrange multiplier - mu, pi0, pi1
u_t1 = -lambda_t1m(4:6);
lambda_t1p46 = lambda_t1p(4:6);
mu_t1p = 1 / (2 * rho^2) * ((r_t1p' * lambda_t1p46) ...
            - (v_t1p' * v_t1p) - (r_t1p' * M1 * r_t1p) ...
            - (r_t1p' * M2 * v_t1p));

pi0 = (-(2 * mu_t1p) * (v_t1p' * v_t1p) - (2 * mu_t1p) * (r_t1p' * ...
      (M1 * r_t1p + M2 + v_t1p + u_t1))) / (4 * (r_t1p' * v_t1p));
pi1 = mu_t1p;

% Equality constraint
res_vec = [2 * (pi0 * r_t1p + pi1 * v_t1p); 2 * pi1 * r_t1p];
ceq = norm(lambda_t1p - lambda_t1m - res_vec);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ----------------------- ODE - dynamics equations -----------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
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
mu = 1 / (2 * rho^2) * (r' * lambda46 - v' * v - r' * M1 * r - r' * M2 * v);

dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v + 2 * mu * r - lambda46;
dydt(7:9) = (4 * mu) * M1 * r - M1 * lambda46 - (2 * mu) * M2 * v + (4 * mu^2) * r - (2 * mu) * lambda46;
dydt(10:12) = M2 * lambda46 - lambda13 + (4 * mu) * v - (2 * mu) * M2 * r;
end