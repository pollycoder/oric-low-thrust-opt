%-------------------------------------------------------------------%
% Spherically constraint - Indirect (Penalty)                       %
% Main Function                                                     %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear;

%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
rho = 10; 
rho0 = 10; omega = 4;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

% Time
t0 = 0;
tf = 0.25;


%-------------------------------------------------------------------%
%------------------------ Guess and Mesh ---------------------------%
%-------------------------------------------------------------------%
n = 50;
tmesh = linspace(t0, tf, n);
yguess = ones(12, 1);
solinit = bvpinit(tmesh, yguess);
options = bvpset('Stats','on','RelTol',1e-5, 'Nmax', 100000);

tic
sol = bvp4c(@bvpfun_eq_penalty, @bvpbc_eq, solinit, options);
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

% Integral
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

save data\indirect_eq_penalty_data.mat rho0 rho x y z ...
                               u1 u2 u3 r u ...
                               tSolve lambda t J
