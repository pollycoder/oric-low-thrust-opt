%-------------------------------------------------------------------%
% Indirect method  - Interior point constraint  (1 point)           %
% LEO: omega = 4 rad/h                                              %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear
tol = 1e-7;
X0 = [0.1248; 0.1252; 1; 1; 1; 1; -8504.99354753256; 
    26299.0055978846; 529.522575635948;	-1863.45274011531;	
    1119.54130372363;	72.6163140659369];
J0 = inf;
while(true)
    [X, J] = fminsearch(@obj_func, X0); 
    if J < J0
        J0 = J;
        X0 = X;
        fprintf("J1=%f\n",J0);
    end
    if J0 < tol
        break
    end
end

J0 = obj_func_en(X0);

n = 0;
while(true)
    [X, J] = fminsearch(@obj_func_en, X0); 
    res = obj_func(X);
    if J < J0
        J0 = J;
        X0 = X;
        fprintf("J2=%f\n",J0); 
    end
    n = n + 1;
    if n > 200
        break
    end
end



%%
tol = 1e-7; J = 0;

%--------------------------- Constant ------------------------------%
t0 = 0;
tf = 0.25;
omega = 4;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

%-------------------- Initial and final state ----------------------%
theta0 = pi;
thetaf = theta0 + pi/2;
rho0 = 10;
rho = 8;
x0 = [rho0*cos(theta0); rho0*sin(theta0); 0; 0; 0; pi];
xf = [rho0*cos(thetaf); rho0*sin(thetaf); 0; 0; 0; pi];

%------------------------ Interior points --------------------------%
t1 = X(1);
scale = 0.1;

%------------------------- Costate guess ---------------------------%
lambda0 = zeros(6, 1);
lambda0(:) = X(7:12);

%--------------------------- Time span -----------------------------%
tspan01 = [t0, t1];
tspan12 = [t1, t2];
tspan2f = [t2, tf];

%----------------------------- Solve -------------------------------%
% Arc1: t0 ~ t1, off
y0 = [x0; lambda0];
[tlist01, y01] = ode45(@odefun_off, tspan01, y0);
x1 = y01(end, 1:6)';
lambda1_m = y01(end, 7:12)';
u01 = -y01(:, 10:12)';
I = trapz(tlist01, 0.5.* (vecnorm(u01).^2));
J =J + I;
res1 = abs(sum(x1(1:3).^2) - rho^2);      % B.C. 1: r(t1) = 8 
res2 =  abs(dot(x1(1:3), x1(4:6)));       % B.C. 2: r(t1)*v(t1) = 0 
    

% Arc2: t1 ~ tf, off
lambda1_p = zeros(6, 1);
lambda1_p(1:3) = lambda1_m(1:3) + (2*pi0) .* x1(1:3) ...
                                + (2*pi1) .* x1(4:6);
lambda1_p(4:6) = lambda1_m(4:6) + (2*pi1) .* x1(1:3);
y1p = [x1, lambda1_p];
[tlist12, y12] = ode45(@odefun_on, tspan12, y1p);
r = y12(:, 1:3)';
v = y12(:, 4:6)';
lambda46 = y12(:, 10:12)';
for i=1:size(r, 2)
    mu(i) = 1 / (2 * rho^2) * (r(:, i)' * lambda46(:, i) ...
                - v(:, i)' * v(:, i) - r(:, i)' * M1 * r(:, i) ...
                - r(:, i)' * M2 * v(:, i));
end
for i=1:size(r, 2)
    u12(:, i) = 2 * mu(i) * r(:, i) - lambda46(:, i);
end
I = trapz(tlist12, 0.5.* (vecnorm(u12).^2));
J = J + I;

r = [vecnorm(y01(:, 1:3)'), vecnorm(y12(:, 1:3)')];
t = [tlist01; tlist12; tlist2f];
plot(t, r);




%-------------------------------------------------------------------%
%---------------------- Objective Function -------------------------%
%-------------------------------------------------------------------%
function J = obj_func(X)
% X = [t1, t2, pi0, pi1, pi2, pi3, l1~l6]
% "*_m" means "-", "*_p" means "+"
J = 0;

%--------------------------- Constant ------------------------------%
t0 = 0;
tf = 0.25;
omega = 4;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

%-------------------- Initial and final state ----------------------%
theta0 = pi;
thetaf = theta0 + pi/2;
rho0 = 10;
rho = 8;
x0 = [rho0*cos(theta0); rho0*sin(theta0); 0; 0; 0; pi];
xf = [rho0*cos(thetaf); rho0*sin(thetaf); 0; 0; 0; pi];

%------------------------ Interior points --------------------------%
t1 = X(1); t2 = X(2);
scale = 0.1;
pi0 = X(3)*scale; pi1 = X(4)*scale; 
pi2 = X(5)*scale; pi3 = X(6)*scale;

%------------------------- Costate guess ---------------------------%
lambda0 = zeros(6, 1);
lambda0(:) = X(7:12);

%--------------------------- Time span -----------------------------%
tspan01 = [t0, t1];
tspan12 = [t1, t2];
tspan2f = [t2, tf];

%----------------------------- Solve -------------------------------%
% Arc1: t0 ~ t1, off
y0 = [x0; lambda0];
[~, y01] = ode45(@odefun_off, tspan01, y0);
x1 = y01(end, 1:6)';
lambda1_m = y01(end, 7:12)';

res1 = abs(sum(x1(1:3).^2) - rho^2);      % B.C. 1: r(t1) = 8 
res2 =  abs(dot(x1(1:3), x1(4:6)));       % B.C. 2: r(t1)*v(t1) = 0 
    
% Arc2: t1 ~ t2, on
lambda1_p = zeros(6, 1);
lambda1_p(1:3) = lambda1_m(1:3) + (2*pi0) .* x1(1:3) ...
                                + (2*pi1) .* x1(4:6);
lambda1_p(4:6) = lambda1_m(4:6) + (2*pi1) .* x1(1:3);
y1p = [x1, lambda1_p];
[~, y12] = ode45(@odefun_on, tspan12, y1p);

% Arc3: t2 ~ tf, off
x2 = y12(end, 1:6)';
lambda2_m = y12(end, 7:12)';
lambda2_p = zeros(6, 1);
lambda2_p(1:3) = lambda2_m(1:3) +  (2*pi2) .* x2(1:3) ...
                                + (2*pi3) .* x2(4:6);
lambda2_p(4:6) = lambda1_m(4:6) + (2*pi3) .* x2(1:3);
y2 = [x2; lambda2_p];
[~, y2f] = ode45(@odefun_off, tspan2f, y2);
xf_solve = y2f(end, 1:6)';
res3 = norm(xf_solve - xf);

J = res1 + res2 + res3;

if t2 < t1
    J = 1e20;
end

end

function J = obj_func_en(X)
% X = [t1, t2, pi0, pi1, pi2, pi3, l1~l6]
% "*_m" means "-", "*_p" means "+"
tol = 1e-7; J = 0;

%--------------------------- Constant ------------------------------%
t0 = 0;
tf = 0.25;
omega = 4;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

%-------------------- Initial and final state ----------------------%
theta0 = pi;
thetaf = theta0 + pi/2;
rho0 = 10;
rho = 8;
x0 = [rho0*cos(theta0); rho0*sin(theta0); 0; 0; 0; pi];
xf = [rho0*cos(thetaf); rho0*sin(thetaf); 0; 0; 0; pi];

%------------------------ Interior points --------------------------%
t1 = X(1); t2 = X(2);
scale = 0.1;
pi0 = X(3)*scale; pi1 = X(4)*scale; 
pi2 = X(5)*scale; pi3 = X(6)*scale;

%------------------------- Costate guess ---------------------------%
lambda0 = zeros(6, 1);
lambda0(:) = X(7:12);

%--------------------------- Time span -----------------------------%
tspan01 = [t0, t1];
tspan12 = [t1, t2];
tspan2f = [t2, tf];

%----------------------------- Solve -------------------------------%
% Arc1: t0 ~ t1, off
y0 = [x0; lambda0];
[tlist01, y01] = ode45(@odefun_off, tspan01, y0);
x1 = y01(end, 1:6)';
lambda1_m = y01(end, 7:12)';
u01 = -y01(:, 10:12)';
I = trapz(tlist01, 0.5.* (vecnorm(u01).^2));
J =J + I;
res1 = abs(sum(x1(1:3).^2) - rho^2);      % B.C. 1: r(t1) = 8 
res2 =  abs(dot(x1(1:3), x1(4:6)));       % B.C. 2: r(t1)*v(t1) = 0 
    

% Arc2: t1 ~ t2, on
lambda1_p = zeros(6, 1);
lambda1_p(1:3) = lambda1_m(1:3) + (2*pi0) .* x1(1:3) ...
                                + (2*pi1) .* x1(4:6);
lambda1_p(4:6) = lambda1_m(4:6) + (2*pi1) .* x1(1:3);
y1p = [x1, lambda1_p];
[tlist12, y12] = ode45(@odefun_on, tspan12, y1p);
r = y12(:, 1:3)';
v = y12(:, 4:6)';
lambda46 = y12(:, 10:12)';
for i=1:size(r, 2)
    mu(i) = 1 / (2 * rho^2) * (r(:, i)' * lambda46(:, i) ...
                - v(:, i)' * v(:, i) - r(:, i)' * M1 * r(:, i) ...
                - r(:, i)' * M2 * v(:, i));
end
for i=1:size(r, 2)
    u12(:, i) = 2 * mu(i) * r(:, i) - lambda46(:, i);
end
I = trapz(tlist12, 0.5.* (vecnorm(u12).^2));
J = J + I;

% Arc3: t2 ~ tf, off
x2 = y12(end, 1:6)';
lambda2_m = y12(end, 7:12)';
lambda2_p = zeros(6, 1);
lambda2_p(1:3) = lambda2_m(1:3) +  (2*pi2) .* x2(1:3) ...
                                + (2*pi3) .* x2(4:6);
lambda2_p(4:6) = lambda1_m(4:6) + (2*pi3) .* x2(1:3);
y2 = [x2; lambda2_p];
[tlist2f, y2f] = ode45(@odefun_off, tspan2f, y2);
xf_solve = y2f(end, 1:6)';
u2f = -y2f(:, 10:12)';
I = trapz(tlist2f, 0.5.* (vecnorm(u2f).^2));
J = J + I;
res3 = norm(xf_solve - xf);

if res1 + res2 + res3 > tol || t2 < t1
    J = 1e20;
end
end











%-------------------------------------------------------------------%
%---------------------- Dynamics Equationn -------------------------%
%-------------------------------------------------------------------%

%----------------------- Unconstrained arc -------------------------%
function dydt = odefun_off(t, y)
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

%------------------------ Constrained arc --------------------------%
function dydt = odefun_on(t, y)
dydt = zeros(12, 1);

% Constant
omega = 4;                                  % angular velocity, 4 rad/h
rho = 8;

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




