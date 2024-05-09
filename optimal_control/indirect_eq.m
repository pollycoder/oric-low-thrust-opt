%-------------------------------------------------------------------%
% Indirect method - Equality Constraint                             %
% LEO: omega = 4 rad/h                                              %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear

%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
omega = 4;                                   % angular velocity, 4 rad/h
rho0 = 10;                                   % Distance between chief and deputy
rho = 10;

% Initial and final states
theta0 = pi;
thetaf = theta0 + pi / 2;
x0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];
xf = [rho0 * cos(thetaf); rho0 * sin(thetaf); 0; 0; 0; pi];

% Time
t0 = 0;
tf = 0.25;

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
A = [zeros(3), eye(3); M1, M2];
B = [zeros(3); eye(3)];


%-------------------------------------------------------------------%
%------------------------ Guess and Mesh ---------------------------%
%-------------------------------------------------------------------%
n = 1e2;
tmesh = linspace(t0, tf, n);
yguess = ones(12, 1);
solinit = bvpinit(tmesh, yguess);
options = bvpset('Stats','on','RelTol',1e-12, 'Nmax', 100000);

tic
sol = bvp5c(@bvpfun_eq, @bvpbc_eq, solinit, options);
tSolve = toc;


%-------------------------------------------------------------------%
%------------------ Calculation and Save Data ----------------------%
%-------------------------------------------------------------------%
r = sol.y(1:3, :);
v = sol.y(4:6, :);
lambda13 = sol.y(7:9, :);
lambda46 = sol.y(10:12, :);
lambda = sol.y(7:12, :);
mu = zeros(size(r, 2), 1);
tol = 1e-5;

% Lagrange multiplier
for i=1:size(r, 2)
    C = norm(r);
    mu(i) = 1 / (2 * rho^2) * (r(:, i)' * lambda46(:, i) ...
                - v(:, i)' * v(:, i) - r(:, i)' * M1 * r(:, i) ...
                - r(:, i)' * M2 * v(:, i));
end

% Integration
for i=1:size(r, 2)
    u(:, i) = 2 * mu(i) * r(:, i) - lambda46(:, i);
end
dJ = 0.5 .* (vecnorm(u) .* vecnorm(u));
t = sol.x;
J = trapz(t, dJ);
fprintf('J = %f', J);

% State - radius
x = sol.y(1, :);
y = sol.y(2, :);
z = sol.y(3, :);
r = sqrt(x.^2 + y.^2 + z.^2);

% Control
u1 = u(1,:);
u2 = u(2,:);
u3 = u(3,:);
u = sqrt(u1.^2 + u2.^2 + u3.^2);

save data\indirect_eq_data.mat rho0 rho x y z u1 u2 u3 r u mu tSolve lambda t J

