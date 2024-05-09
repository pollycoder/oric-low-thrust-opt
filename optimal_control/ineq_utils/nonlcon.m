%-------------------------------------------------------------------%
% Indirect method  - Interior point constraint  (2 points)          %
% Objective Function                                                %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
function [c, ceq] = nonlcon(X)
%-------------------------------------------------------------------%
% "*_p" <==> "t+"                                                   %
% "*_m" <==> "t-"                                                   %
% X: dt1, dt2, lambda0(6), lambda1+(6), lambda2+(6),                %
%    theta1, phi1, theta2, phi2, v1(3), v2(3)                       %
%-------------------------------------------------------------------%

%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
rho = 9; 
rho0 = 10;
theta0 = pi; thetaf = theta0 + pi/2;
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
[~, y01] = ode45(@odefun_off, [t0, t1], y0_guess);
y1m = y01(end, :)';

% Arc1: on
[~, y12] = ode45(@odefun_on, [t1, t2], y1p_guess);
y2m = y12(end, :)';
y1p = y12(1, :)';

% Arc3: off
[~, y2f] = ode45(@odefun_off, [t2, tf], y2p_guess);
y2p = y2f(1, :)';
yf = y2f(end, :)';

%---------------------------- Residuals -----------------------------%
% Res1: final state
res1 = yf(1:6) - statef;

% Res2: joint points jump and state continuity - t1
res2 = [y1p(1:6)-y1m(1:6); y1p(10:12)-y1m(10:12)];

% Res3: joint points jump and state continuity - t2
res3 = [y2p(1:6)-y2m(1:6); y2p(10:12)-y2m(10:12)];

% Res4: tangent condition
res4 = [dot(r1guess, v1guess); dot(r2guess, v2guess)];

% Constriant
x = [y01(:,1);y12(:,1);y2f(:,1)];
y = [y01(:,2);y12(:,2);y2f(:,2)];
z = [y01(:,3);y12(:,3);y2f(:,3)];
r = sqrt(x.^2 + y.^2 + z.^2);
rbase = rho .* ones(size(r));

% Final residual
ceq = [res1; res2; res3; res4];

c = max(rbase - r);

end

