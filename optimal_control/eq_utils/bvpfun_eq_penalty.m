%-------------------------------------------------------------------%
% Spherically constraint - Indirect (Pontryagin)                    %
% Dynamics Equation                                                 %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
function dydt = bvpfun_eq_penalty(t, y)
% Constant
omega = 4;                                  % angular velocity, 4 rad/h
alpha = 1e6;                                % Parameter to be adjusted
rho = 10;                                   % Distance between chief and deputy

% Matrix
M1 = diag([3 * omega^2, 0, -omega^2]);
M2 = diag([2 * omega, 0], 1) + diag([-2 * omega, 0], -1);
A = [zeros(3), eye(3); M1, M2];
B = [zeros(3); eye(3)];

syms x1 x2 x3 p1 p2 p3
C =x1^2 + x2^2 + x3^2 - rho^2;
C_value = double(subs(C, {x1, x2, x3}, {y(1), y(2), y(3)}));
gradC = [diff(C, x1); diff(C, x2); diff(C, x3); diff(C, p1); diff(C, p2); diff(C, p3)];
gradC_value = double(subs(gradC, {x1, x2, x3, p1, p2, p3}, {y(1), y(2), y(3), y(4), y(5), y(6)}));
dydt(1:6) = A * y(1:6) - B * B' * y(7:12);
dydt(7:12) = -A' * y(7:12) - alpha * C_value * gradC_value;
end





