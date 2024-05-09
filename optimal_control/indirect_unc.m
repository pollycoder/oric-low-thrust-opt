%-------------------------------------------------------------------%
% Indirect method - Unconstrained Problem                           %
% Main function                                                     %
% LEO: omega = 4 rad/h                                              %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear

tic

% Constant
omega = 4;                                  % angular velocity, 4 rad/h
rho = 10;                                   % Distance between chief and deputy

% Initial and final states
x0 = [-rho; 0; 0; 0; 0; pi];
xf = [0; -rho; 0; 0; 0; pi];
t0 = 0;
tf = 0.25   ;

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
yguess = 10*ones(12, 1);
solinit = bvpinit(tmesh, yguess);
options = bvpset('Stats','on','RelTol',1e-9, 'Nmax', 100000);

sol = bvp4c(@bvpfun_unc, @bvpbc_unc, solinit, options);

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
title('Trajectory');


