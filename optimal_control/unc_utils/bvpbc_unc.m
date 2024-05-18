%-------------------------------------------------------------------%
% Unconstrained - Indirect (Pontryagin)                             %
% Boundary Conditions                                               %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
function res = bvpbc_unc(y0, yf)
rho0 = 10;
theta0 = pi;
thetaf = theta0 + pi / 2;
x0 = [rho0 * cos(theta0); rho0 * sin(theta0); 9; 8; 2; pi];
xf = [rho0 * cos(thetaf); rho0 * sin(thetaf); 3; 1; 3; pi];

xf = 2.*x0;

% Initial and final states
res = [y0(1:6) - x0; yf(1:6) - xf];
end