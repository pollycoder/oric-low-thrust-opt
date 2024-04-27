%-------------------------------------------------------------------%
% Spherically constrained - GPOPS-II                                %
% Endpoint function                                                 %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
function output = sphereOptEndpoint(input)
q = input.phase.integral;
output.objective = q;

end

