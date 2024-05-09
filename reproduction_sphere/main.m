%-------------------------------------------------------------------%
% Spherically constrained - CubicFitRot Reproduction                %
% Main function                                                     %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
clc;clear; tic

%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
syms t
x0 = [-10, 0, 0, 0, 0, pi]';
xf = [0, -10, 0, 0, 0, pi]';
t0 = 0;
tf = 0.25;
tValue = linspace(t0, tf, 1000);
g = 100;
omega = 4;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
s0 = Cartesian2Spherical(x0);
sf = Cartesian2Spherical(xf);


%-------------------------------------------------------------------%
%-------------------------- CubicFitRot ----------------------------%
%-------------------------------------------------------------------%
[T, u, s, J, rho] = cubicFitRot(x0, xf, t0, tf, M1, M2, g);
tSolve = toc;

u = double(subs(u, t, tValue));
s4D = double(subs(s, t, tValue'));

n = size(s4D, 1);
s6D = zeros(n, 6);
s6D(:, 1) = rho * ones(1, n);
s6D(:, 4) = zeros(1, n);
s6D(:, [2, 3, 5, 6]) = s4D;

[x, y, z] = sph2cart(s6D(:, 2), pi / 2 * ones(n, 1) - s6D(:, 3), s6D(:, 1));
X = [x';y';z'];
XBack = T \ X;
x = XBack(1, :)';
y = XBack(2, :)';
z = XBack(3, :)';


%-------------------------------------------------------------------%
%------------------------------- Plot ------------------------------%
%-------------------------------------------------------------------%
figure
plot3(x, y, z, 'LineWidth', 2, 'Color', 'k'); hold on
xlabel('x');
ylabel('y');
zlabel('z');

r = rho;
[X, Y, Z] = sphere;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2, Y2, Z2,  'FaceAlpha', 0.2, 'EdgeColor', ...
    'texturemap'); hold on
colormap(gca, 'bone')
axis equal

plot3(0, 0, 0, '*k', 'LineWidth', 2);hold on
text(0, 0, 0, 'Chief'); hold on
plot3(x(1), y(1), z(1), 'vg', 'LineWidth', 2); hold on
text(x(1), y(1), z(1), 'Deputy-Departure');hold on
plot3(x(end), y(end), z(end), 'vb', 'LineWidth', 2);hold on
text(x(end), y(end), z(end), 'Deputy-Arrival');hold on

index = linspace(1, length(x), 112);
x1 = x(index);
y1 = y(index);
z1 = z(index);
u1x = u(1, index)';
u1y = u(2, index)';
u1z = u(3, index)';

figure
u = u1x.^2 + u1y.^2 + u1z.^2;
plot(tValue(index), u, 'LineWidth', 1.5);




