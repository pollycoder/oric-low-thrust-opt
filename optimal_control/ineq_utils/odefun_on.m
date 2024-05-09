%-------------------------------------------------------------------%
% Indirect method  - Interior point constraint  (2 points)          %
% Motion equations - "on" arc                                       %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
function dydt = odefun_on(t, y)
dydt = zeros(12, 1);

% Constant
omega = 4;                                  % angular velocity, 4 rad/h
rho = 9;

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

% Lagrange multiplier
r = y(1:3);
v = y(4:6);
lambda13 = y(7:9);
lambda46 = y(10:12);

eta = 1/(2*rho^2) * (2*v'*M2*lambda46 - 4*v'*lambda13 ...
                   - 3*(lambda46'*lambda46) + 8*r'*M1*lambda46 ...
                   + 2*r'*M2*M2*lambda46 - 2*r'*M2*lambda13 ...
                   -4*v'*M1*v - 4*r'*M1*M1*r - 3*r'*M1*M2*v ...
                   - r'*M2*M1*v);

dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v - lambda46;
dydt(7:9) = - M1 * lambda46 + (2 * eta) * r;
dydt(10:12) = M2 * lambda46 - lambda13;

end


