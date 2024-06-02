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
rho = 9.5; 
rho0 = 10; omega = 4;
theta0 = pi;
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

t0 = 0; tf = 0.25;
dt1 = X(1); dt2 = X(2);
t1 = t0 + dt1; t2 = t1 + dt2;
state0 = [rho0 * cos(theta0); rho0 * sin(theta0); 0; 0; 0; pi];


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

%-------------------------- Integration ----------------------------%
t = [t01; t12; t2f];
x = [y01(:,1); y12(:,1); y2f(:,1)];
y = [y01(:,2); y12(:,2); y2f(:,2)];
z = [y01(:,3); y12(:,3); y2f(:,3)];
r = sqrt(x.^2 + y.^2 + z.^2);

% Costate
costate = [y01(:,7:12)', y12(:,7:12)', y2f(:,7:12)'];
lambda = costate;

% Control
uvec = -[y01(:, 10:12)', y12(:, 10:12)', y2f(:, 10:12)'];

% Lagrange multiplier - mu and eta
lambda46 = lambda(4:6, :);
state = [y01(:,1:6)', y12(:,1:6)', y2f(:,1:6)'];
mu = zeros(length(r), 1);
for i=1:length(r)
    if i <= size(y01, 1) || i > size(y01, 1) + size(y12, 1)
        mu(i) = 0;
    else
        mu(i) = 1 / (2 * rho^2) * (state(1:3, i)' * lambda46(:, i) ...
                - state(4:6, i)' * state(4:6, i) ...
                - state(1:3, i)' * M1 * state(1:3, i) ...
                - state(1:3, i)' * M2 * state(4:6, i));
    end
    uvec(:, i) = 2 .* mu(i) .* [x(i); y(i); z(i)] + uvec(:, i);
end

% Control norm
u = vecnorm(uvec);

% Integration target
dJ = 0.5 * u.^2;
J = trapz(t, dJ);

Jhat = J;

end




