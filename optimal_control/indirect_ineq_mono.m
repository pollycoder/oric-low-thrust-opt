%-------------------------------------------------------------------%
% Indirect method  - Interior point constraint  (1 point)           %
% LEO: omega = 4 rad/h                                              %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear









%-------------------------------------------------------------------%
%---------------------- Objective Function -------------------------%
%-------------------------------------------------------------------%
function obj_func(X)

% X = [t1, theta1]
t0 = 0;
tf = 0.25;
omega = 4;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

% Initial and final state
theta0 = pi;
thetaf = theta0 + pi/2;
rho0 = 10;
rho = 8;
x0 = [rho0*cos(theta0); rho*sin(theta0); 0; 0; 0; pi];
xf = [rho0*cos(thetaf); rho*sin(thetaf); 0; 0; 0; pi];

% Interior points (n=2)
r1 = [rho*cos(theta1); rho*sin(theta1); 0];


end


%-------------------------------------------------------------------%
%---------------------- Dynamics Equationn -------------------------%
%-------------------------------------------------------------------%

%----------------------- Unconstrained arc -------------------------%
function dydt = bvpfun_off(t, y)
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
mu = 0;

% Equations
dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v + 2 * mu * r - lambda46;
dydt(7:9) = (4 * mu) * M1 * r - M1 * lambda46 ...
            - (2 * mu) * M2 * v + (4 * mu^2) * r - (2 * mu) * lambda46;
dydt(10:12) = M2 * lambda46 - lambda13 + (4 * mu) * v - (2 * mu) * M2 * r;
end


%-------------------------------------------------------------------%
%------------------------------- Bounds ----------------------------%
%-------------------------------------------------------------------%
function res = bvpbc(y0, yf)
% Initial and final states
rho0 = 10;
x0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];
xf = [rho0 * cos(thetaf); rho0 * sin(thetaf); 0; 0; 0; pi];
res = [y0(1:6) - x0; yf(1:6) - xf];
end



