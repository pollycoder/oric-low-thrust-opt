clc
clear
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Indirect method
% LEO: omega = 4 rad/h
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global rho rho0 x0 xf
% Constant
omega = 4;                                   % angular velocity, 4 rad/h
rho0 = 10;                                   % Distance between chief and deputy
rho = 11.79868;

% Initial and final states
theta0 = pi;
thetaf = theta0 + pi / 2;
x0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];
xf = [rho0 * cos(thetaf); rho0 * sin(thetaf); 0; 0; 0; pi];

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
n = 10;
tmesh = linspace(t0, tf, n);
yguess = ones(12, 1);
solinit = bvpinit(tmesh, yguess);
options = bvpset('Stats','on','RelTol',1e-9, 'Nmax', 100000);

tic
sol = bvp5c(@bvpfun, @bvpbc, solinit, options);
tSolve = toc;

%%
r = sol.y(1:3, :);
v = sol.y(4:6, :);
lambda13 = sol.y(7:9, :);
lambda46 = sol.y(10:12, :);
lambda = sol.y(7:12, :);
mu = zeros(size(r, 2), 1);
tol = 1e-5;

for i=1:size(r, 2)
    C = norm(r);
    mu(i) = 1 / (2 * rho^2) * (r(:, i)' * lambda46(:, i) ...
                - v(:, i)' * v(:, i) - r(:, i)' * M1 * r(:, i) ...
                - r(:, i)' * M2 * v(:, i));
end


for i=1:size(r, 2)
    u(:, i) = 2 * mu(i) * r(:, i) - lambda46(:, i);
end
dJ = 0.5 .* (vecnorm(u) .* vecnorm(u));
t = sol.x;
J = trapz(t, dJ);
fprintf('J = %f', J);


%% State - radius
x = sol.y(1, :);
y = sol.y(2, :);
z = sol.y(3, :);
r = sqrt(x.^2 + y.^2 + z.^2);

%% Plot
figure
set(gca, 'XTick', 9:0.5:11);
plot(t, r, 'LineWidth', 1.5);
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
u1 = u(1,:);
u2 = u(2,:);
u3 = u(3,:);
u = sqrt(u1.^2 + u2.^2 + u3.^2);
plot(t, u1, 'LineWidth', 1.5);hold on
plot(t, u2, 'LineWidth', 1.5);hold on
plot(t, u3, 'LineWidth', 1.5);
legend('control1', 'control2', 'control3');
title('Control');

figure
plot(t, u, 'LineWidth', 1.5);
title('Control - Norm');

% Multiplier
figure
plot(t, mu, 'LineWidth', 1.5);
title('mu');

% Trajectory
figure
rb = rho;
index = 1:1000:size(u, 2);
[X, Y, Z] = sphere;
X2 = X * rb;
Y2 = Y * rb;
Z2 = Z * rb;
surf(X2, Y2, Z2,  'FaceAlpha', 0.2, 'EdgeColor', 'texturemap'); hold on
colormap(gca, 'bone')
axis equal
plot3(0, 0, 0, 'k*', 'LineWidth', 3);hold on
text(0, 0, 0, 'Chief');hold on
plot3(x(1), y(1), z(1), 'g*', 'LineWidth', 2);hold on
text(x(1), y(1), z(1), 'Departure');hold on
plot3(x(end), y(end), z(end), 'c*', 'LineWidth', 2);hold on
text(x(end), y(end), z(end), 'Arrival');hold on
plot3(x, y, z, 'k-', 'LineWidth', 1.5);hold on
title('Trajectory');

save data\indirect_lag_data.mat x y z u1 u2 u3 r u mu tSolve

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
mu = 1 / (2 * rho^2) * (r' * lambda46 - v' * v - r' * M1 * r - r' * M2 * v);

dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v + 2 * mu * r - lambda46;
dydt(7:9) = (4 * mu) * M1 * r - M1 * lambda46 - (2 * mu) * M2 * v + (4 * mu^2) * r - (2 * mu) * lambda46;
dydt(10:12) = M2 * lambda46 - lambda13 + (4 * mu) * v - (2 * mu) * M2 * r;
end

% Boundary conditions
function res = bvpbc(y0, yf)
global rho0 x0 xf
% Initial and final states
res = [y0(1:6) - x0; yf(1:6) - xf];
end