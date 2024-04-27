%-------------------------------------------------------------------%
% Spherically constrained - GPOPS-II                                %
% Main function                                                     %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
clc;clear

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
rho_ub = 15;

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
vmax = 10;
vmin = -vmax;
umax = 2000;
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

%------------------------------ Guess ------------------------------%
% Time
guess.phase.time = [t0; tf];

% State
guess.phase.state = [x0; x0];

% Control
guess.phase.control(:, 1) = [770; 780];
guess.phase.control(:, 2) = [1050; 1100];
guess.phase.control(:, 3) = [-100; -50];

% Integral
guess.phase.integral = 2e5;


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
setup.mesh.method = 'hp1';
setup.mesh.tolerance = 1e-6; 
setup.mesh.phase.colpoints = 4*ones(1,10);
setup.mesh.phase.fraction = 0.1*ones(1,10);


%-------------------------------------------------------------------%
%------------------- Solve Problem Using GPOPS2 --------------------%
%-------------------------------------------------------------------%
output = gpops2(setup);
output.result.nlptime
solution = output.result.solution;
