%-------------------------------------------------------------------%
% Spherically Inequality Constraint - Indirect (Pontryagin)         %
% Main Function                                                     %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear
%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
rho0 = 10; rho = 9;
theta0 = pi; omega = 4;
thetaf = theta0 + pi / 2;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
x0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];

% Time
t0 = 0;
tf = 0.25;


%-------------------------------------------------------------------%
%------------------------ Guess and Mesh ---------------------------%
%-------------------------------------------------------------------%
n = 1e2;
tmesh = linspace(t0, tf, n);
yguess = [x0; 1e3.*ones(6, 1)];
solinit = bvpinit(tmesh, yguess);
options = bvpset('Stats','on','RelTol',1e-12, 'Nmax', 100000);

tic
sol = bvp5c(@bvpfun_sf, @bvpbc_sf, solinit, options);
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
                             - v(:, i)' * v(:, i) ...
                             - r(:, i)' * M1 * r(:, i) ...
                             - r(:, i)' * M2 * v(:, i));
end

% Integral
for i=1:size(r, 2)
    u(:, i) = -lambda46(:, i);
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

% Lagrange multiplier - mu and eta
lambda13 = lambda(1:3, :);
lambda46 = lambda(4:6, :);
state = sol.y(1:6, :);
mu = zeros(length(r), 1);
eta1 = zeros(length(r), 1);
eta2 = zeros(length(r), 1);
for i=1:length(r)
    eta1(i) = 1/(2*rho^2) * (2*state(4:6, i)'*M2*lambda46(:, i) ...
                   - 4*state(4:6, i)'*lambda13(:, i) ...
                   - 3*(lambda46(:, i)'*lambda46(:, i)) ...
                   + 8*state(1:3, i)'*M1*lambda46(:, i) ...
                   + 2*state(1:3, i)'*M2*M2*lambda46(:, i) ...
                   - 2*state(1:3, i)'*M2*lambda13(:, i) ...
                   -4*state(1:3, i)'*M1*state(1:3, i) ...
                   - 4*state(1:3, i)'*M1*M1*state(1:3, i) ...
                   - 3*state(1:3, i)'*M1*M2*state(4:6, i) ...
                   - state(1:3, i)'*M2*M1*state(4:6, i));
    eta2(i) = 1/(2*rho^2) * (3*state(4:6, i)'*lambda46(:, i) ...
                           - state(1:3, i)'*lambda13(:, i) ...
                           + state(1:3, i)'*M2*lambda46(:, i) ...
                           -4*state(4:6, i)'*M1*state(1:3, i) ...
                           - state(1:3, i)'*M2*M1*state(1:3, i) ...
                           - state(1:3, i)'*M2*M2*state(4:6, i) ...
                           + state(1:3, i)'*M2*lambda46(:, i));
    mu(i) = 1 / (2 * rho^2) * (state(1:3, i)' * lambda46(:, i) ...
                - state(4:6, i)' * state(4:6, i) ...
                - state(1:3, i)' * M1 * state(1:3, i) ...
                - state(1:3, i)' * M2 * state(4:6, i));
end
    
save data/indirect_ineq_sf_data.mat rho0 rho x y z ...
                               u1 u2 u3 r u mu ...
                               tSolve lambda t J mu eta1 eta2