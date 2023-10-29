clc
clear workspace
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Reproduction of spherically-constrained optimal control
% Paper resource:
%   Woodford, N. T., Harris, M. W., & Petersen, C. D. (2023). 
% Spherically constrained relative motion trajectories in 
% low earth orbit. Journal of Guidance, Control, and Dynamics, 
% 46(4), 666-679.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
syms t

x0 = [-10, 0, 0, 0, 0, pi]';
xf = [0, -10, 0, 0, 0, pi]';
t0 = 0;
tf = 0.25;
tValue = linspace(t0, tf, 1000);

omega = 4;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
[T, u, s, J, rho] = cubicFitRot(x0, xf, t0, tf, M1, M2, 3);
TBack = blkdiag(inv(T), inv(T));


u = double(subs(u, t, tValue));
s4D = double(subs(s, t, tValue'));

n = size(s4D, 1);
s6D = zeros(n, 6);
s6D(:, 1) = rho * ones(1, n);
s6D(:, 4) = zeros(1, n);
s6D(:, [2, 3, 5, 6]) = s4D;

[x, y, z] = sph2cart(s6D(:, 2), pi / 2 * ones(n, 1) - s6D(:, 3), s6D(:, 1));


plot3(x, y, z, 'LineWidth', 2, 'Color', 'b'); hold on
xlabel('x');
ylabel('y');
zlabel('z');

r = 8;
[X, Y, Z] = sphere;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2, Y2, Z2,  'FaceAlpha', 0.5, 'EdgeColor', 'texturemap'); hold on
colormap(gca, 'copper')
axis equal

plot3(0, 0, 0, '*k', 'LineWidth', 2);hold on
plot3(x(1), y(1), z(1), '*g', 'LineWidth', 2); hold on
text(x(1), y(1), z(1), 'Departure');hold on
plot3(x(end), y(end), z(end), '*r', 'LineWidth', 2);hold on
text(x(end), y(end), z(end), 'Arrival');hold on

index = linspace(1, length(x), 112);
x1 = x(index);
y1 = y(index);
z1 = z(index);
u1x = u(1, index)';
u1y = u(2, index)';
u1z = u(3, index)';
quiver3(x1, y1, z1, u1x, u1y, u1z, 0.3, 'Color', 'r','LineWidth', 1);







