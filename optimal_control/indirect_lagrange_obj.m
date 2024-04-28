%-------------------------------------------------------------------%
% Indirect method                                                   %
% LEO: omega = 4 rad/h                                              %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear

%-------------------------------------------------------------------%
%--------------------------- Optimization --------------------------%
%-------------------------------------------------------------------%
rho_init = 11.78;
[rho, J] = fminsearch(@obj_func, rho_init);


%-------------------------------------------------------------------%
%------------------------ Objective Function -----------------------%
%-------------------------------------------------------------------%
function J = obj_func(rho)

%---------------------------- Constant -----------------------------%
omega = 4;                                   % angular velocity, 4 rad/h
rho_min = 8;

% Time
t0 = 0;
tf = 0.25;

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

%------------------------- Guess and Mesh --------------------------%
n = 1e2;
tmesh = linspace(t0, tf, n);
yguess = ones(12, 1);
solinit = bvpinit(tmesh, yguess);
options = bvpset('Stats','on','RelTol',1e-12, 'Nmax', 100000);

%-------------------- Calculation and Save Data --------------------%
sol = bvp4c(@(t, y)bvpfun(t, y, rho), @bvpbc, solinit, options);

r = sol.y(1:3, :);
v = sol.y(4:6, :);
lambda46 = sol.y(10:12, :);
mu = zeros(size(r, 2), 1);

% Lagrange multiplier
for i=1:size(r, 2)
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

rNorm = vecnorm(r);
rmin = min(rNorm);

if rmin < rho_min
    J = 1e6;
end

end


%-------------------------------------------------------------------%
%---------------------- Dynamics Equationn -------------------------%
%-------------------------------------------------------------------%
function dydt = bvpfun(t, y, rho)
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
mu = 1 / (2 * rho^2) * (r' * lambda46 - v' * v ...
                        - r' * M1 * r - r' * M2 * v);

% Equations
dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v + 2 * mu * r - lambda46;
dydt(7:9) = (4 * mu) * M1 * r - M1 * lambda46 ...
            - (2 * mu) * M2 * v + (4 * mu^2) * r - (2 * mu) * lambda46;
dydt(10:12) = M2 * lambda46 - lambda13 + (4 * mu) * v - (2 * mu) * M2 * r;
end


%-------------------------------------------------------------------%
%------------------------------ Bounds -----------------------------%
%-------------------------------------------------------------------%
function res = bvpbc(y0, yf)
rho0 = 10;
theta0 = pi;
thetaf = theta0 + pi / 2;
x0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];
xf = [rho0 * cos(thetaf); rho0 * sin(thetaf); 0; 0; 0; pi];

% Initial and final states
res = [y0(1:6) - x0; yf(1:6) - xf];
end