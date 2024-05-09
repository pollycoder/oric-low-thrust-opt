%-------------------------------------------------------------------%
% Spherically constrained - GPOPS-II                                %
% Dynamics equation                                                 %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
function phaseout = sphereOptContinuous(input)
M1 = input.auxdata.M1;
M2 = input.auxdata.M2;

t = input.phase.time;
x = input.phase.state;
u = input.phase.control;

dx(:, 1:3) = x(:, 4:6);
dx(:, 4:6) = x(:, 1:3) * M1' + x(:, 4:6) * M2' + u;

phaseout.dynamics = dx;
phaseout.integrand = 1/2 * (u(:,1).^2 + u(:,2).^2 + u(:,3).^2);
phaseout.path = sqrt(x(:, 1).^2 + x(:, 2).^2 + x(:, 3).^2);

end

