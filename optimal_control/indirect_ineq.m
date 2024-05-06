%-------------------------------------------------------------------%
% Indirect method  - Interior point constraint  (1 point)           %
% LEO: omega = 4 rad/h                                              %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear
global rho rho0 x0 xf

%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
omega = 4;                                   % angular velocity, 4 rad/h
rho0 = 10;                                   % Distance between chief and deputy
rho = 8;

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
sol = bvp5c(@bvpfun, @bvpbc, solinit, options);
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
    mu(i) = 0;
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

save data\indirect_ineq_data.mat rho0 rho x y z u1 u2 u3 r u mu tSolve lambda t J


%-------------------------------------------------------------------%
%---------------------- Dynamics Equationn -------------------------%
%-------------------------------------------------------------------%
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
eta = (1/(2*rho^2)) * (r' * lambda13 + r' * M2 * M2 * lambda46 ...
                       + r' * M1 * lambda46 - r' * M2 * lambda13);

% Equations
dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v - lambda46;
dydt(7:9) = - M1 * lambda46 + (2 * eta) * r;
dydt(10:12) = M2 * lambda46 - lambda13;
end


%-------------------------------------------------------------------%
%------------------------------ Bounds -----------------------------%
%-------------------------------------------------------------------%
function res = bvpbc(y0, yf)
global rho0 x0 xf
% Initial and final states
res = [y0(1:6) - x0; yf(1:6) - xf];
end