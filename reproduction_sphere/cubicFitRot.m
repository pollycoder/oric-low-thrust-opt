%-------------------------------------------------------------------%
% Spherically constrained - CubicFitRot Reproduction                %
% CubicFitRot Algorithm                                             %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
function [TStar, uStar, sStar, JStar, rho] = cubicFitRot(CartesianState0, ...
                                                         CartesianStatef, ...
                                                         t0, tf, M1, M2, g, tol)
%-------------------------------------------------------------------%
% Input:                                                            %
%   s0: Initial state - theta0, phi0, d(theta)/dt(0), d(phi)/dt(0)  %
%   sf: Final state - thetaf, phif, d(theta)/dt(f), d(phi)/dt(f)    %
%   t0: Initial time                                                %
%   tf: Final time                                                  %
%   M1,M2: C-W matrix                                               %
%   rho: radial length                                              %
%   g - the number of search angles                                 %
%-------------------------------------------------------------------%
% Output:                                                           %
%   J: Optimizing index                                             %
%   u: controlled variable                                          %
%   state: theta, phi, d(theta)/dt, d(phi)/dt                       %
%   T - the rotating matrix                                         %
%-------------------------------------------------------------------%

if nargin < 8
    tol = 1e-3;
end
if nargin < 7
    g = 3;
end
if nargin < 6
    error("Not enough inputs. See cubicFitRot.");
end
ang = linspace(0, 2*pi, g);
JStar = inf;
TStar = zeros(6);
uStar = zeros(3, 1);
sStar = zeros(6, 1);

for i=1:g
    for j=1:g
        % Matrix T
        alpha1 = ang(i);
        alpha2 = ang(j);
        T1 = diag([1, cos(alpha1), cos(alpha1)]) + ...
             diag([0, sin(alpha1)], -1) + diag([0, -sin(alpha1)], 1);
        T2 = diag([cos(alpha2), 1, cos(alpha2)]) + ...
             diag(-sin(alpha2), -2) + diag(sin(alpha2), 2);
        T = T1 * T2;
        
        % Rotate State
        RVHat0 = blkdiag(T, T) * CartesianState0;
        RVHatf = blkdiag(T, T) * CartesianStatef;

        % M1 M2
        M1Hat = (T * M1) / T;
        M2Hat = (T * M2) / T;

        % Coordinate transform
        s0Hat6D = Cartesian2Spherical(RVHat0);
        sfHat6D = Cartesian2Spherical(RVHatf);
        
        % Variable extracting
        err = s0Hat6D(1) - sfHat6D(1);
        temp = abs(err);
        if temp > tol
            error("Radial coordinate not the same, can't satisfy spherical " + ...
                "constraint.");
        end
        rho = s0Hat6D(1);
        s0Hat4D = s0Hat6D([2, 3, 5, 6]);
        sfHat4D = sfHat6D([2, 3, 5, 6]);

        % Cubic fit
        [J, uHat, sHat4D] = cubicFit(s0Hat4D, sfHat4D, t0, tf, M1Hat, M2Hat, rho);
        if J < JStar
            JStar = J;
            TStar = T;
            uStar = T \ uHat;
            sStar = sHat4D;
        end        
    end
end

end