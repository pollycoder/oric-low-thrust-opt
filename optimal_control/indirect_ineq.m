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
rho = 9; 
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
lambda0 = [22.9444423724889	
           33530.8874064454	
           463.372437424655	
           -1486.67431726342	
           1613.82304610368	
           69.1580010840004];
lambda1 = [-16450.1815123402	
            19001.2973265684	
            673.740161518210	
           -1414.68332810443	
           -717.708453657806
            18.5528906459845];
lambda2 = [-30850.6541415830	
           -3700.71005509594	
            616.707806830944	
           -859.126418227609	
           -1002.92293539350	
           -14.7458926735603];

% t1 and t2
t1 = 0.0998167723302513;
t2 = 0.149074508752072;

% Velocities at t1 and t2
v1 = [41.6148287058002; -64.7193778424496; -1.36570613231753];
v2 = [62.6476324259827; -42.7456628770326; -1.56291874984878];

% Positions at t1 and t2 (spherical coordinate)
theta1 = -2.55949814018984;
theta2 = -2.14145326036684;
phi1 = 0.00505190120742112;
phi2 = -0.00364951776919052;

% Initial values
X0 = [t1; t2-t1; lambda0;lambda1;lambda2; 
      theta1; phi1; theta2; phi2; v1; v2];

%------------------- Parameters for optimization -------------------%
tol = 0.04;
A = zeros(3, 30);
A(1,1) = -1;
A(2,2) = -1;
A(3, 1:2) = ones(1, 2);
b = zeros(3, 1);
b(1) = -tol;
b(2) = -tol;
b(3) = tf - t0;
options = optimoptions("fmincon", ...
                       "ConstraintTolerance", 1e-12, ...
                       "FunctionTolerance", 1e-12, ...
                       "MaxIterations", 1e6, ...
                       "UseParallel",true, ...
                       "MaxFunctionEvaluations", 5e5, ...
                       "StepTolerance", 1e-15);
[X, J] = fmincon(@obj_func, X0, A, b, [],[],[],[], @nonlcon,options);


%-------------------------------------------------------------------%
%------------------- Resolve the original problem ------------------%
%-------------------------------------------------------------------%

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
v1guess = X(index6:index7-1);
v2guess = X(index7:index8-1);

% State guess
state1_guess = [r1guess; v1guess];
state2_guess = [r2guess; v2guess];

%------------------------ Multiple Shooting ------------------------%
tic
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

%--------------------------- Residuals -----------------------------%
% Res1: final state
res1 = yf(1:6) - statef;

% Res2: joint points jump and state continuity - t1
res2 = [y1p(1:6)-y1m(1:6); y1p(10:12)-y1m(10:12)];

% Res3: joint points jump and state continuity - t2
res3 = [y2p(1:6)-y2m(1:6); y2p(10:12)-y2m(10:12)];

% Res4: tangent condition
res4 = [dot(r1guess, v1guess); dot(r2guess, v2guess)];


%-------------------------------------------------------------------%
%--------------------------- Save Data -----------------------------%
%-------------------------------------------------------------------%
t = [t01; t12; t2f];
x = [y01(:,1); y12(:,1); y2f(:,1)];
y = [y01(:,2); y12(:,2); y2f(:,2)];
z = [y01(:,3); y12(:,3); y2f(:,3)];
r = sqrt(x.^2 + y.^2 + z.^2);

% Costate
costate = [y01(:,7:12)', y12(:,7:12)', y2f(:,7:12)'];
lambda = costate;

% Control
u01 = -y01(:, 10:12)';
u12 = -y12(:, 10:12)';
u2f = -y2f(:, 10:12)';
u1 = [u01(1,:), u12(1,:), u2f(1,:)]';
u2 = [u01(2,:), u12(2,:), u2f(2,:)]';
u3 = [u01(3,:), u12(3,:), u2f(3,:)]';
u01norm = vecnorm(u01);
u12norm = vecnorm(u12);
u2fnorm = vecnorm(u2f);
u = [u01norm, u12norm, u2fnorm]';

% Integration target
dJ01 = 0.5*u01norm.^2;
dJ12 = 0.5*u12norm.^2;
dJ2f = 0.5*u2fnorm.^2;
J = trapz(t01, dJ01) + trapz(t12, dJ12) + trapz(t2f, dJ2f);

% Lagrange multiplier - mu and eta
lambda13 = lambda(1:3, :);
lambda46 = lambda(4:6, :);
state = [y01(:,1:6)', y12(:,1:6)', y2f(:,1:6)'];
mu = zeros(length(r), 1);
eta = zeros(length(r), 1);
for i=1:length(r)
    if i > size(y01, 1) && i <= size(y01,1) + size(y12,1)
        eta(i) = 1/(2*rho^2) * (2*state(4:6, i)'*M2*lambda46(:, i) ...
                   - 4*state(4:6, i)'*lambda13(:, i) ...
                   - 3*(lambda46(:, i)'*lambda46(:, i)) ...
                   + 8*state(1:3, i)'*M1*lambda46(:, i) ...
                   + 2*state(1:3, i)'*M2*M2*lambda46(:, i) ...
                   - 2*state(1:3, i)'*M2*lambda13(:, i) ...
                   -4*state(1:3, i)'*M1*state(1:3, i) ...
                   - 4*state(1:3, i)'*M1*M1*state(1:3, i) ...
                   - 3*state(1:3, i)'*M1*M2*state(4:6, i) ...
                   - state(1:3, i)'*M2*M1*state(4:6, i));
    else
        eta(i) = 0;
    end
    mu(i) = 1 / (2 * rho^2) * (state(1:3, i)' * lambda46(:, i) ...
                - state(4:6, i)' * state(4:6, i) ...
                - state(1:3, i)' * M1 * state(1:3, i) ...
                - state(1:3, i)' * M2 * state(4:6, i));
end

save data\indirect_ineq_data.mat rho0 rho x y z ...
                                 u1 u2 u3 r u tSolve ...
                                 lambda t J mu eta






















