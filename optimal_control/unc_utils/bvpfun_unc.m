%-------------------------------------------------------------------%
% Indirect method - Unconstrained Problem                           %
% Boundary Conditions                                               %
% LEO: omega = 4 rad/h                                              %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
function dydt = bvpfun_unc(t, y)
dydt = zeros(12, 1);
% Constant
omega = 4;                                  % angular velocity, 4 rad/h

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);

r = y(1:3);
v = y(4:6);
lambda13 = y(7:9);
lambda46 = y(10:12);

dydt(1:3) = v;
dydt(4:6) = M1 * r + M2 * v - lambda46;
dydt(7:9) = - M1 * lambda46;
dydt(10:12) = M2 * lambda46 - lambda13;
end