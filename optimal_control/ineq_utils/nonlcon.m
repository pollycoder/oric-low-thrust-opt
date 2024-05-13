%-------------------------------------------------------------------%
% Spherically Inequality Constraint - Indirect (Pontryagin)         %
% Non-linear Constraint for fmincon                                 %
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
rho0 = 10; omega = 4;
theta0 = pi; thetaf = theta0 + pi/2;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
t0 = 0; tf = 0.25;
state0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];
statef = [rho0 * cos(thetaf); rho0 * sin(thetaf); 0; 0; 0; pi];

dt1 = X(1); dt2 = X(2);
t1 = t0 + dt1; t2 = t1 + dt2;


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

% Data
x = [y01(:,1); y12(:,1); y2f(:,1)];
y = [y01(:,2); y12(:,2); y2f(:,2)];
z = [y01(:,3); y12(:,3); y2f(:,3)];
r = sqrt(x.^2 + y.^2 + z.^2);

% Costate
costate = [y01(:,7:12)', y12(:,7:12)', y2f(:,7:12)'];
lambda = costate;

% Lagrange multiplier - mu and eta
lambda13 = lambda(1:3, :);
lambda46 = lambda(4:6, :);
state = [y01(:,1:6)', y12(:,1:6)', y2f(:,1:6)'];
mu = zeros(length(r), 1);
eta1 = zeros(length(r), 1);
eta2 = zeros(length(r), 1);
for i=1:length(r)
    if i > size(y01, 1) && i <= size(y01,1) + size(y12,1)
        
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
    else
        eta1(i) = 0;
        eta2(i) = 0;
    end
    mu(i) = 1 / (2 * rho^2) * (state(1:3, i)' * lambda46(:, i) ...
                - state(4:6, i)' * state(4:6, i) ...
                - state(1:3, i)' * M1 * state(1:3, i) ...
                - state(1:3, i)' * M2 * state(4:6, i));
end
%mu1_p = mu(size(y01, 1)+1);

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
ceq = [res1; res2; res3];%; res4];
c = [max(rbase - r)];

end

