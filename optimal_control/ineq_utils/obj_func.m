%-------------------------------------------------------------------%
% Spherically Inequality Constraint - Indirect (Pontryagin)         %
% Objective Function                                                %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
function Jhat = obj_func(X)
%-------------------------------------------------------------------%
% "*_p" <==> "t+"                                                   %
% "*_m" <==> "t-"                                                   %
% X: dt1, dt2, lambda0(6), lambda1+(6), lambda2+(6),                %
%    theta1, phi1, theta2, phi2, v1(3), v2(3)                       %
%-------------------------------------------------------------------%
X = reshape(X, [30, 1]);
%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%
rho = 9; 
rho0 = 10;
theta0 = pi;

t0 = 0; tf = 0.25;
dt1 = X(1); dt2 = X(2);
t1 = t0 + dt1; t2 = t1 + dt2;
state0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];


%-------------------------------------------------------------------%
%--------------------------- Constraint ----------------------------%
%-------------------------------------------------------------------%
alpha_eq = 1e4;
alpha_neq = 1e3;
ctol = 0;
%----------------------------- Linear ------------------------------%
tol = 0.04;
A = zeros(3, 30);
A(1,1) = -1;
A(2,2) = -1;
A(3, 1:2) = ones(1, 2);
b = zeros(3, 1);
b(1) = -tol;
b(2) = -tol;
b(3) = tf - t0;
res1 = sign(max(A*X - b) + ctol) + 1;


%--------------------------- Non Linear ----------------------------%
[c, ceq] = nonlcon(X);
res2 = sum(ceq.^2);
res3 = sign(max(c) + ctol) + 1;

penalty = alpha_eq * res2 + alpha_neq * (res1 + res3);


%-------------------------------------------------------------------%
%-------------------------- Costate Guess --------------------------%
%-------------------------------------------------------------------%
index0 = 3; dimLambda = 6; 
index1 = index0 + dimLambda;
index2 = index1 + dimLambda;
index3 = index2 + dimLambda;

lambda0 = zeros(dimLambda, 1);
lambda1p_guess = zeros(dimLambda, 1);
lambda2p_guess = zeros(dimLambda, 1);

for i=1:dimLambda
    lambda0(i) = X(index0+i-1);
    lambda1p_guess(i) = X(index1+i-1);
    lambda2p_guess(i) = X(index2+i-1);
end


%-------------------------------------------------------------------%
%--------------------------- State Guess ---------------------------%
%-------------------------------------------------------------------%
dimSphere = 2; dimV = 3;
index4 = index3;
index5 = index4 + dimSphere;
index6 = index5 + dimSphere;
index7 = index6 + dimV;

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


%-------------------------------------------------------------------%
%--------------------------- State Guess ---------------------------%
%-------------------------------------------------------------------%
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-10, ...
                 'NormControl','on');

%------------------------ Multiple Shooting ------------------------%
y0_guess = [state0; lambda0];
y1p_guess = [state1_guess; lambda1p_guess];
y2p_guess = [state2_guess; lambda2p_guess];

% Arc1: off
[t01, y01] = ode45(@odefun_off, [t0, t1], y0_guess, options);

% Arc1: on
[t12, y12] = ode45(@odefun_on, [t1, t2], y1p_guess, options);

% Arc3: off
[t2f, y2f] = ode45(@odefun_off, [t2, tf], y2p_guess, options);

%------------------------------ Index ------------------------------%
u01 = -y01(:, 10:12)';
u12 = -y12(:, 10:12)';
u2f = -y2f(:, 10:12)';
u01norm = vecnorm(u01);
u12norm = vecnorm(u12);
u2fnorm = vecnorm(u2f);
dJ01 = 0.5*u01norm.^2;
dJ12 = 0.5*u12norm.^2;
dJ2f = 0.5*u2fnorm.^2;

%penalty = 0;
Jhat = trapz(t01, dJ01) + trapz(t12, dJ12) + trapz(t2f, dJ2f)...
       + penalty ;

end




