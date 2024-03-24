clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Indirect method
% LEO: omega = 4 rad/h
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic
global rho
% Constant
omega = 4;                                  % angular velocity, 4 rad/h
%rho0 = 10;
rho = 10;                                   % Distance between chief and deputy

% Initial and final states
x0 = [-rho; 0; 0; 0; 0; pi];
xf = [0; -rho; 0; 0; 0; pi];
t0 = 0;
tf = 0.25;

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
A = [zeros(3), eye(3); M1, M2];
B = [zeros(3); eye(3)];

% Initial guess for the solution 2  
% t - [t0, tf]
% x - x0, [1,1,1,1,1,1]
n = 1e3;
tmesh = linspace(t0, tf, n);
yguess = ones(12, 1);
solinit = bvpinit(tmesh, yguess);
options = bvpset('Stats','on','RelTol',1e-9, 'Nmax', 100000);

sol = bvp4c(@bvpfun, @bvpbc, solinit, options);

r = sol.y(1:3, :);
v = sol.y(4:6, :);
lambda13 = sol.y(7:9, :);
lambda46 = sol.y(10:12, :);
lambda = sol.y(7:12, :);
mu = zeros(size(r, 2), 1);
for i=1:size(r, 2)
    C = norm(r);
    mu(i) = 0;
end

for i=1:size(r, 2)
    u(:, i) = 2 * mu(i) * r(:, i) - lambda46(:, i);
end
dJ = 0.5 .* (vecnorm(u) .* vecnorm(u));
t = sol.x;
J = trapz(t, dJ);
fprintf('J = %f', J);
tSolve = toc;

%% State - radius
x1 = sol.y(1, :);
x2 = sol.y(2, :);
x3 = sol.y(3, :);
x = sqrt(x1.^2 + x2.^2 + x3.^2);
x = x - rho * ones(size(x));

%% Plot
figure
set(gca, 'XTick', 9:0.5:11);
plot(t, x, 'LineWidth', 1.5);
title('State - pos');

% Costate
costate = lambda;
figure
plot(t, costate(1, :), 'LineWidth', 1.5);hold on
plot(t, costate(2, :), 'LineWidth', 1.5);hold on
plot(t, costate(3, :), 'LineWidth', 1.5);hold on
plot(t, costate(4, :), 'LineWidth', 1.5);hold on
plot(t, costate(5, :), 'LineWidth', 1.5);hold on
plot(t, costate(6, :), 'LineWidth', 1.5);
legend('costate1', 'costate2', 'costate3', 'costate4', 'costate5', 'costate6');
title('Costate');

% Control
figure
plot(t, sqrt(2 .* dJ), 'LineWidth', 1.5);
title('Control');

% Trajectory
figure
r = rho;
index = 1:1000:size(u, 2);
uIndex = u(:, index);
x1Index = x1(index);
x2Index = x2(index);
x3Index = x3(index);
[X, Y, Z] = sphere;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2, Y2, Z2,  'FaceAlpha', 0.2, 'EdgeColor', 'texturemap'); hold on
colormap(gca, 'bone')
axis equal
plot3(0, 0, 0, 'k*', 'LineWidth', 3);hold on
text(0, 0, 0, 'Chief');hold on
plot3(x1(1), x2(1), x3(1), 'g*', 'LineWidth', 2);hold on
text(x1(1), x2(1), x3(1), 'Departure');hold on
plot3(x1(end), x2(end), x3(end), 'c*', 'LineWidth', 2);hold on
text(x1(end), x2(end), x3(end), 'Arrival');hold on
plot3(x1, x2, x3, 'k-', 'LineWidth', 1.5);hold on
%quiver3(x1Index, x2Index, x3Index, uIndex(1, :), uIndex(2, :), uIndex(3, :), 0.3, 'Color', 'r','LineWidth', 1.5);
title('Trajectory');
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Control equations
% f = dydt
% y = [x; lambda], 12x1 vector
% x: state - v and a, 6x1 vector
% lambda: costate, 6x1 vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function dydt = bvpfun(t, y)
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

% Boundary conditions
function res = bvpbc(y0, yf)
global rho
% Initial and final states
x0 = [-rho; 0; 0; 0; 0; pi];
xf = [0; -rho; 0; 0; 0; pi];
res = [y0(1:6) - x0; yf(1:6) - xf];
end