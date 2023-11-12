clc
clear all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Indirect method
% LEO: omega = 4 rad/h
% P3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global omega alpha rho A B x0 xf
% Constant
omega = 4;                                  % angular velocity, 4 rad/h
alpha = 1e5;                                 % Parameter to be adjusted
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
tmesh = linspace(t0, tf, 50);
yguess = [8.521; -4.894; -1.854; -5.529; -8.154; -3.385; ones(6, 1)];
solinit = bvpinit(tmesh, yguess);
options = bvpset('Stats','on','RelTol',1e-5);

sol = bvp4c(@bvpfun, @bvpbc, solinit, options);

% Plot
r = 10;
[X, Y, Z] = sphere;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2, Y2, Z2,  'FaceAlpha', 0.2, 'EdgeColor', 'texturemap'); hold on
colormap(gca, 'bone')
axis equal

x1 = sol.y(1, :);
x2 = sol.y(2, :);
x3 = sol.y(3, :);
x = x1.^2 + x2.^2 + x3.^2;

figure
plot3(0, 0, 0, 'k*', 'LineWidth', 3);hold on
plot3(x1(1), x2(1), x3(1), 'g*', 'LineWidth', 2);hold on
plot3(x1(end), x2(end), x3(end), 'r*', 'LineWidth', 2);hold on
plot3(x1, x2, x3, 'k-', 'LineWidth', 1.5);

% Control equations
% f = dydt
% y = [x; lambda], 12x1 vector
% x: state - v and a, 6x1 vector
% lambda: costate, 6x1 vector
function dydt = bvpfun(t, y)
global rho alpha A B
syms x1 x2 x3 p1 p2 p3
C =x1^2 + x2^2 + x3^2 - rho^2;
C_value = double(subs(C, {x1, x2, x3}, {y(1), y(2), y(3)}));
gradC = [diff(C, x1); diff(C, x2); diff(C, x3); diff(C, p1); diff(C, p2); diff(C, p3)];
gradC_value = double(subs(gradC, {x1, x2, x3, p1, p2, p3}, {y(1), y(2), y(3), y(4), y(5), y(6)}));
dydt(1:6) = A * y(1:6) - B * B' * y(7:12);
dydt(7:12) = -A' * y(7:12) - alpha * C_value * gradC_value;
end

% Boundary conditions
function res = bvpbc(y0, yf)
global x0 xf
res = [y0(1:6) - x0; yf(1:6) - xf];
end




