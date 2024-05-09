%-------------------------------------------------------------------%
% Indirect method - Equality Constraint                             %
% Dynamics Equation                                                 %
% LEO: omega = 4 rad/h                                              %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
function dydt = bvpfun_eq(t, y)
rho = 10;
dydt = zeros(12, 1);

% Constant
omega = 4;                                  % angular velocity, 4 rad/h

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

% Lagrange multiplier
r = y(1:3);
v = y(4:6);
lambda13 = y(7:9);
lambda46 = y(10:12);
mu = 1 / (2 * rho^2) * (r' * lambda46 - v' * v ...
                        - r' * M1 * r - r' * M2 * v);

% Equations
dydt(1:3) = v;
dydt(4:6) = M1 * y(1:3) + M2 * v + 2 * mu * r - lambda46;
dydt(7:9) = (4 * mu) * M1 * r - M1 * lambda46 ...
            - (2 * mu) * M2 * v + (4 * mu^2) * r - (2 * mu) * lambda46;
dydt(10:12) = M2 * lambda46 - lambda13 + (4 * mu) * v - (2 * mu) * M2 * r;
end


