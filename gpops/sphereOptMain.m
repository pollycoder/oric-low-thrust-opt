%-------------------------------------------------------------------%
% Spherically constrained - GPOPS-II                                %
% Main function                                                     %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
clc;clear
tic
%-------------------------------------------------------------------%
%---------------------------- Constant -----------------------------%
%-------------------------------------------------------------------%

% Matrix for dynamics equation
omega = 4;
auxdata.M1 = diag([3 * omega^2, 0, -omega^2]);
auxdata.M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

% Distance - initial, final, bounds
rho0 = 10;
rhof = 10;
rho_lb = 8;
rho_ub = 10;

% Time - initial, final
t0 = 0;
tf = 0.25;

% State - initial, final
theta0 = pi;
thetaf = theta0 + pi/2;
v0 = [0, 0, pi];
vf = [0, 0, pi];
x0 = [rho0*[cos(theta0), sin(theta0), 0], v0];
xf = [rho0*[cos(thetaf), sin(thetaf), 0], vf];

% Bounds
rmax = rho_ub;
rmin = -rmax;
vmax = 65;
vmin = -vmax;
umax = 1400;
umin = -umax;


%-------------------------------------------------------------------%
%------------------------ Bounds and Guess -------------------------%
%-------------------------------------------------------------------%

%----------------------------- Bounds ------------------------------%
% Time
bounds.phase.initialtime.lower = t0;
bounds.phase.initialtime.upper = t0;
bounds.phase.finaltime.lower = tf;
bounds.phase.finaltime.upper = tf;

% State
bounds.phase.initialstate.lower = x0;
bounds.phase.initialstate.upper = x0;
bounds.phase.finalstate.lower = xf;
bounds.phase.finalstate.upper = xf;
bounds.phase.state.lower = [rmin*ones(1,3), vmin*ones(1,3)];
bounds.phase.state.upper = [rmax*ones(1,3), vmax*ones(1,3)];

% Control
bounds.phase.control.lower = umin * ones(1, 3);
bounds.phase.control.upper = umax * ones(1, 3);

% Path constraint
bounds.phase.lower = rho_lb;
bounds.phase.upper = rho_ub;

% Integral
bounds.phase.integral.lower = 0;
bounds.phase.integral.upper = 5e5;
bounds.phase.path.lower = rho_lb;
bounds.phase.path.upper = rho_ub;

%------------------------------ Guess ------------------------------%
% Use spherical coordinate to guess position
N = 100;
thetaGuess = zeros(N, 1);
thetaGuess(:) = linspace(theta0, thetaf, N);
rhoGuess = zeros(N, 1);
tol = (rho0 - rho_lb) * 2 / N;
rhoGuess(1:N/2) = linspace(rho0, rho_lb, N/2);
rhoGuess(N/2+1:N) = linspace(rho_lb + tol, rhof, N/2);

% Velocity guess
vmiddle = [vmax, vmin, -2];
vtol = (vmiddle - v0) * 2 ./ N;
vGuess = zeros(N, 3);
vGuess(1:N/2, 1) = linspace(v0(1), vmiddle(1), N/2);
vGuess(N/2+1:N, 1) = linspace(vmiddle(1)+vtol(1), vf(1), N/2);
vGuess(1:N/2, 2) = linspace(v0(2), vmiddle(2), N/2);
vGuess(N/2+1:N, 1) = linspace(vmiddle(2)+vtol(2), vf(2), N/2);
vGuess(1:N/2, 3) = linspace(v0(3), vmiddle(3), N/2);
vGuess(N/2+1:N, 3) = linspace(vmiddle(3)+vtol(3), vf(3), N/2);

% Time
guess.phase.time = linspace(t0, tf, N)';

% State
guess.phase.state(:, 1) = rhoGuess .* cos(thetaGuess);
guess.phase.state(:, 2) = rhoGuess .* sin(thetaGuess);
guess.phase.state(:, 3) = rhoGuess .* zeros(N, 1);
guess.phase.state(:, 4) = vGuess(:, 1);
guess.phase.state(:, 5) = vGuess(:, 2);
guess.phase.state(:, 6) = vGuess(:, 3);

% Control
guess.phase.control(:, 1) = linspace(umax, umin, N);
guess.phase.control(:, 2) = linspace(umin, umax, N);
guess.phase.control(:, 3) = linspace(-70, 70, N);

% Integral
guess.phase.integral = 0;


%-------------------------------------------------------------------%
%-------------------------------- Mesh -----------------------------%
%-------------------------------------------------------------------%
mesh.method = 'hp1';
mesh.tolerance = 1e-9; 
mesh.maxiteration = 1000;
mesh.colpointmin = 4;
mesh.colpointmax = 50;
mesh.phase.colpoints = 4*ones(1,50);
mesh.phase.fraction = 0.02*ones(1,50);



%-------------------------------------------------------------------%
%--------------------------- Problem Setup -------------------------%
%-------------------------------------------------------------------%
setup.name = 'Spherically Constrained Optimization';
setup.functions.continuous = @sphereOptContinuous;
setup.functions.endpoint = @sphereOptEndpoint;
setup.auxdata = auxdata;
setup.bounds = bounds;
setup.guess = guess;

setup.nlp.solver = 'snopt';
setup.derivatives.supplier = 'sparseCD';
setup.derivatives.derivativelevel = 'second';
setup.scales.method = 'automatic-bounds';
setup.mesh = mesh;
setup.method = 'RPMintegration';


%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
tSolve = output.result.nlptime;
solution = output.result.solution;


%%
t = solution.phase.time;

figure
x = solution.phase.state(:, 1);
y = solution.phase.state(:, 2);
z = solution.phase.state(:, 3);
r = sqrt(x.^2 + y.^2 + z.^2);
plot(t, r, 'LineWidth', 1.5);
title('State - pos');

figure
u1 = solution.phase.control(:, 1);
u2 = solution.phase.control(:, 2);
u3 = solution.phase.control(:, 3);
u = sqrt(u1.^2 + u2.^2 + u3.^2);
plot(t, u1, 'LineWidth', 1.5);hold on
plot(t, u2, 'LineWidth', 1.5);hold on
plot(t, u3, 'LineWidth', 1.5);hold on
legend('u1', 'u2', 'u3');
title('Control');

figure
plot(t, u, 'LineWidth', 1.5);
title('Control - Norm');

% Trajectory
figure
rb = rho0;
index = 1:1000:size(u, 2);
[X, Y, Z] = sphere;
X2 = X * rb;
Y2 = Y * rb;
Z2 = Z * rb;
surf(X2, Y2, Z2,  'FaceAlpha', 0.2, 'EdgeColor', 'texturemap'); hold on
colormap(gca, 'bone')
axis equal
plot3(0, 0, 0, 'k*', 'LineWidth', 3);hold on
text(0, 0, 0, 'Chief');hold on
plot3(x(1), y(1), z(1), 'g*', 'LineWidth', 2);hold on
text(x(1), y(1), z(1), 'Departure');hold on
plot3(x(end), y(end), z(end), 'c*', 'LineWidth', 2);hold on
text(x(end), y(end), z(end), 'Arrival');hold on
plot3(x, y, z, 'k-', 'LineWidth', 1.5);hold on
title('Trajectory');

save data\gpops_data.mat x y z u1 u2 u3 r u tSolve