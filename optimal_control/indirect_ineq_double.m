%-------------------------------------------------------------------%
% Indirect method  - Interior point constraint  (1 point)           %
% LEO: omega = 4 rad/h                                              %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear

lambda0 = [-8488.56386370120	
            26299.8437937099	
            529.325246545735	
            -1862.30789185224	
            1119.73400284263	
            72.6080715242726];
lambda1 = [-17140.5777965151	
            10387.3458218917	
            609.744314834945	
            -1143.63936992369	
            -893.419700734222	
            0.612161804832503];
lambda2 = [-27573.7462572614	
            27.4214092241810	
            612.596813195079	
            -1136.30407075007	
            -891.282009812638	
            0.331116284423432];

X0 = [0.1243; 0.1244; lambda0;lambda1;lambda2;
      pi*1.25; 0; pi*1.25; 0; 
      50; -50; -1; 50; -50; -1];
t0 = 0; tf = 0.25;

tol = 1e-7;
A = zeros(3, 30);
A(1,1) = -1;
A(2,2) = -1;
A(3, 1:2) = ones(1, 2);
b = zeros(3, 1);
b(1) = -tol;
b(2) = -tol;
b(3) = tf - t0;
options = optimoptions("fmincon", ...
                       "ConstraintTolerance", 1e-8, ...
                       "FunctionTolerance", 1e-8, ...
                       "MaxIterations", 1e5, ...
                       "UseParallel",true, ...
                       "MaxFunctionEvaluations", 4e5);
[X, J] = fmincon(@obj_res, X0, A, b, [],[],[],[], @nonlcon,options);

%%
%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
rho = 8; 
rho0 = 10; omega = 4;
theta0 = pi; thetaf = theta0 + pi/2;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

t0 = 0; tf = 0.25;
dt1 = X(1); dt2 = X(2);
t1 = t0 + dt1; t2 = t1 + dt2;
state0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];
statef = [rho0 * cos(thetaf); rho0 * sin(thetaf); 0; 0; 0; pi];


%-------------------------------------------------------------------%
%-------------------------- Costate Guess --------------------------%
%-------------------------------------------------------------------%
index0 = 3; dimLambda = 6; 
index1 = index0 + dimLambda;
index2 = index1 + dimLambda;
index3 = index2 + dimLambda;
lambda0 = X(index0:index1-1);
lambda1p_guess = X(index1:index2-1);
lambda2p_guess = X(index2:index3-1);


%-------------------------------------------------------------------%
%--------------------------- State Guess ---------------------------%
%-------------------------------------------------------------------%
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


%-------------------------------------------------------------------%
%--------------------------- State Guess ---------------------------%
%-------------------------------------------------------------------%

%------------------------ Multiple Shooting ------------------------%
y0_guess = [state0; lambda0];
y1p_guess = [state1_guess; lambda1p_guess];
y2p_guess = [state2_guess; lambda2p_guess];

% Arc1: off
[t01, y01] = ode45(@odefun_off, [t0, t1], y0_guess);
y1m = y01(end, :)';
lambda1m = y1m(7:12);

% Arc1: on
[t12, y12] = ode45(@odefun_on, [t1, t2], y1p_guess);
y2m = y12(end, :)';
lambda2m = y2m(7:12);

% Arc3: off
[t2f, y2f] = ode45(@odefun_off, [t2, tf], y2p_guess);
yf = y2f(end, :)';

%---------------------------- Residuals -----------------------------%
% Res1: final state
res1 = yf(1:6) - statef;

% Res2: joint points jump and state continuity - t1
r = y1p_guess(1:3);
v = y1p_guess(4:6);
lambda13m = y1m(7:9);
lambda46m = y1m(10:12);
lambda13p = y1p_guess(7:9);
lambda46p = y1p_guess(10:12);
mu1p = 1 / (2 * rho^2) * (r' * lambda46p - v' * v ...
                        - r' * M1 * r - r' * M2 * v);
res2 = [y1p_guess(1:6) - y1m(1:6); 
        lambda13p - lambda13m - mu1p * r;
        lambda46p - lambda46m];

% Res3: joint points jump and state continuity - t2
r = y2p_guess(1:3);
v = y2p_guess(4:6);
lambda13m = y2m(7:9);
lambda46m = y2m(10:12);
lambda13p = y2p_guess(7:9);
lambda46p = y2p_guess(10:12);
mu2p = 1 / (2 * rho^2) * (r' * lambda46m - v' * v ...
                        - r' * M1 * r - r' * M2 * v);
res3 = [y2p_guess(1:6) - y2m(1:6); 
        lambda13p + mu2p * r - lambda13m;
        lambda46p - lambda46m];

% Final residual
res = [res1; res2; res3];



%%
t = [t01; t12; t2f];
x = [y01(:,1);y12(:,1);y2f(:,1)];
y = [y01(:,2);y12(:,2);y2f(:,2)];
z = [y01(:,3);y12(:,3);y2f(:,3)];
r = sqrt(x.^2 + y.^2 + z.^2);
u = -[y01(:, 10:12)', y12(:, 10:12)', y2f(:, 10:12)'];
costate = [y01(:, 7:12)', y12(:, 7:12)', y2f(:, 7:12)'];
unorm = vecnorm(u);
dJ = 0.5*unorm.^2;
J = trapz(t, dJ);

figure
plot(t, r);
title('State')

figure
plot(t, costate(1,:), 'LineWidth', 1.5);hold on
plot(t, costate(2,:), 'LineWidth', 1.5);hold on
plot(t, costate(3,:), 'LineWidth', 1.5);hold on
plot(t, costate(4,:), 'LineWidth', 1.5);hold on
plot(t, costate(5,:), 'LineWidth', 1.5);hold on
plot(t, costate(6,:), 'LineWidth', 1.5);hold on
legend('costate1', 'costate2', 'costate3', 'costate4', ...
       'costate5', 'costate6')
title('Costate')

figure
plot(t, u(1,:), 'LineWidth', 1.5);hold on
plot(t, u(2,:), 'LineWidth', 1.5);hold on
plot(t, u(3,:), 'LineWidth', 1.5);hold on
legend('control1', 'control2', 'control3')
title('Control')

figure
plot(t, unorm);hold on
title('UNorm')






















