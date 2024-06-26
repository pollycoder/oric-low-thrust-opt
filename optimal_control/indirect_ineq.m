%-------------------------------------------------------------------%
% Spherically Inequality Constraint - Indirect (Pontryagin)         %
% Main Function                                                     %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear;tic

%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
rho = 9.5; 
rho0 = 10; omega = 4;
theta0 = pi; thetaf = theta0 + pi/2;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

t0 = 0; tf = 0.25;
state0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];
statef = [rho0 * cos(thetaf); rho0 * sin(thetaf); 0; 0; 0; pi];


%-------------------------------------------------------------------%
%-------------------------- Optimization ---------------------------%
%-------------------------------------------------------------------%

%----------------------------- Guess -------------------------------%
% Lambda at t0, t1, t2
lambda0 = 1e3 .* ones(6, 1);
lambda1 = 1e3 .* ones(6, 1);
lambda2 = 1e3 .* ones(6, 1);

% t1 and t2
t1 = 0.1;
t2 = 0.15;

% Velocities at t1 and t2
v1 = [40 .* ones(2, 1); 1];
v2 = [40 .* ones(2, 1); 1];

% Positions at t1 and t2 (spherical coordinate)
theta1 = -2.6;
theta2 = -2.1;
phi1 = 0.005;
phi2 = -0.004;

% Initial values
X0 = [t1; t2-t1; lambda0;lambda1;lambda2; 
      theta1; phi1; theta2; phi2; v1; v2];

%------------------- Parameters for optimization -------------------%
tol = 1e-6;
A = zeros(3, 30);
A(1,1) = -1;
A(2,2) = -1;
A(3, 1:2) = ones(1, 2);
b = zeros(3, 1);
b(1) = -tol;
b(2) = -tol;
b(3) = tf - t0;

%
options = optimoptions("fmincon", ...
                       "ConstraintTolerance", 1e-8, ...
                       "FunctionTolerance", 1e-8, ...
                       "MaxIterations", 1e3, ...
                       "UseParallel",true, ...
                       "MaxFunctionEvaluations", 1e6, ...
                       "Algorithm", "sqp", ...
                       "OptimalityTolerance", 1e-8, ...
                       "StepTolerance", 1e-15, ...
                       "Display", "iter");
[X, result] = fmincon(@obj_func, X0, A, b, [],[],[],[], @nonlcon, options);
%}

%{
options = optimoptions('fsolve', ...
                       'Algorithm','levenberg-marquardt', ...
                       'FunctionTolerance', 1e-8, ...
                       'OptimalityTolerance', 1e-8, ...
                       'MaxIterations', 1e5, ...
                       'StepTolerance', 1e-15, ...
                       'MaxFunctionEvaluations', 1e5, ...
                       'Display', 'iter', ...
                       'UseParallel',  true);
X = fsolve(@nonlcon, X0, options);
%}

%%
%-------------------------------------------------------------------%
%------------------- Resolve the original problem ------------------%
%-------------------------------------------------------------------%
X = reshape(X, [30, 1]);
%------------------------------- Time ------------------------------%
dt1 = X(1); dt2 = X(2);
t1 = t0 + dt1; t2 = t1 + dt2;

%------------------------------ Costate ----------------------------%
index0 = 3; dimLambda = 6; 
index1 = index0 + dimLambda;
index2 = index1 + dimLambda;
index3 = index2 + dimLambda;
lambda0 = X(index0:index1-1);
lambda1p_guess = X(index1:index2-1);
lambda2p_guess = X(index2:index3-1);

%------------------------------- State -----------------------------%
dimSphere = 2; dimV = 3;
index4 = index3;
index5 = index4 + dimSphere;
index6 = index5 + dimSphere;
index7 = index6 + dimV;
index8 = index7 + dimV;

% Spherical coordinate - angles guess
theta1 = X(index4);
phi1 = X(index5-1);
theta2 = X(index5);
phi2 = X(index6-1);

% Position guess
[x1guess, y1guess, z1guess] = sph2cart(theta1, phi1, rho);
r1guess = [x1guess; y1guess; z1guess];
[x2guess, y2guess, z2guess] = sph2cart(theta2, phi2, rho);
r2guess = [x2guess; y2guess; z2guess];

% Velocity guess
v1guess = zeros(dimV, 1);
v2guess = zeros(dimV, 1);

for i=1:dimV
    v1guess(i) = X(index6+i-1);
    v2guess(i) = X(index7+i-1);
end

% State guess
state1_guess = [r1guess; v1guess];
state2_guess = [r2guess; v2guess];

%------------------------ Multiple Shooting ------------------------%
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-10, ...
                 'NormControl','on');

% y-guess
y0_guess = [state0; lambda0];
y1p_guess = [state1_guess; lambda1p_guess];
y2p_guess = [state2_guess; lambda2p_guess];

% Arc1: off
[t01, y01] = ode45(@odefun_off, [t0, t1], y0_guess, options);
y1m = y01(end, :)';
lambda1m = y1m(7:12);

% Arc1: onqqq
[t12, y12] = ode45(@odefun_on, [t1, t2], y1p_guess, options);
y2m = y12(end, :)';
y1p = y12(1, :)';
lambda2m = y2m(7:12);

% Arc3: off
[t2f, y2f] = ode45(@odefun_off, [t2, tf], y2p_guess, options);
y2p = y2f(1, :)';
yf = y2f(end, :)';
tSolve = toc;


%-------------------------------------------------------------------%
%--------------------------- Save Data -----------------------------%
%-------------------------------------------------------------------%
t = [t01; t12; t2f];
x = [y01(:,1); y12(:,1); y2f(:,1)];
y = [y01(:,2); y12(:,2); y2f(:,2)];
z = [y01(:,3); y12(:,3); y2f(:,3)];
vx = [y01(:,4); y12(:,4); y2f(:,4)];
vy = [y01(:,5); y12(:,5); y2f(:,5)];
vz = [y01(:,6); y12(:,6); y2f(:,6)];
r = sqrt(x.^2 + y.^2 + z.^2);

% Costate
costate = [y01(:,7:12)', y12(:,7:12)', y2f(:,7:12)'];
lambda = costate;

% Control
uvec = -[y01(:, 10:12)', y12(:, 10:12)', y2f(:, 10:12)'];

% Lagrange multiplier - mu and eta
lambda13 = lambda(1:3, :);
lambda46 = lambda(4:6, :);
state = [y01(:,1:6)', y12(:,1:6)', y2f(:,1:6)'];
mu = zeros(length(r), 1);
H = zeros(length(r), 1);
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
    if i <= size(y01, 1) || i > size(y01, 1) + size(y12, 1)
        mu(i) = 0;
        eta1(i) = 0;
    else
        mu(i) = 1 / (2 * rho^2) * (state(1:3, i)' * lambda46(:, i) ...
                - state(4:6, i)' * state(4:6, i) ...
                - state(1:3, i)' * M1 * state(1:3, i) ...
                - state(1:3, i)' * M2 * state(4:6, i));
    end
    H(i) = -1/2 * (lambda46(:, i)' * lambda46(:, i)) + lambda13(:, i)' * state(4:6, i) ...
           + lambda46(:, i)' * M1 * state(1:3, i) + lambda46(:, i)' * M2 * state(4:6, i);
end

% Control components and norm
u1 = uvec(1, :)';
u2 = uvec(2, :)';
u3 = uvec(3, :)';
u = vecnorm(uvec);

% Integration target
dJ = 0.5*u.^2;
J = trapz(t, dJ);


save data/indirect_ineq_data_9_5.mat rho0 rho x y z ...
                                 u1 u2 u3 r u tSolve ...
                                 lambda t J mu eta1 eta2 H


%-------------------------------------------------------------------%
%--------------------------- Residuals -----------------------------%
%-------------------------------------------------------------------%
% Res1: final state
res1 = yf(1:6) - statef;

% Res2: joint points jump and state continuity - t1
res2 = [y1p(1:6)-y1m(1:6); y1p(10:12)-y1m(10:12)];

% Res3: joint points jump and state continuity - t2
res3 = [y2p(1:6)-y2m(1:6); y2p(10:12)-y2m(10:12)];

% Res4: tangent condition
res4 = [dot(r1guess, v1guess); dot(r2guess, v2guess)];
