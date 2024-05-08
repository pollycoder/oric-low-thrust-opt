%-------------------------------------------------------------------%
% Indirect method  - Interior point constraint  (1 point)           %
% LEO: omega = 4 rad/h                                              %
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %  
%-------------------------------------------------------------------%
clc;clear
tol = 1e-7;
lambda0 = [-8488.56386370120	
           26299.8437937099	
           529.325246545735	
           -1862.30789185224	
           1119.73400284263	
           72.6080715242726];

X0 = [0.125; 1e-3; lambda0; lambda0; lambda0;
      pi*5/4; 0; pi*5/4; 0; 
      40; -70; 0.1; 40; -70; -0.1];
res0 = obj_res(X0);
n = 0;
while(true)
    tol = 1e-8;
    options = optimset("UseParallel", true);
    [X, res] = fminsearch(@obj_res, X0, options);
    if res < res0
        X0 = X;
        res0 = res;
        fprintf("Res = %f\n", res0);
    end
    n = n + 1;
    if n > 1000 || res0 < tol
        break;
    end
end

