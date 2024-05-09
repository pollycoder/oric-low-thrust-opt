%-------------------------------------------------------------------%
% Spherically constrained - CubicFitRot Reproduction                %
% Spherical coordinates ---> Cartesian coordinates                  %
%-------------------------------------------------------------------%
% Reference: Woodford N T, Harris M W, Petersen C D. Spherically    %
% constrained relative motion trajectories in low earth orbit[J].   %
% Journal of Guidance, Control, and Dynamics, 2023, 46(4): 666-679. %                                                
%-------------------------------------------------------------------%
function CartesianState = Spherical2Cartesian(SphericalState6D)
%-------------------------------------------------------------------%
% Cartesian states:                                                 %
%   CartesianState(1): x                                            %
%   CartesianState(2): y                                            %
%   CartesianState(3): z                                            %
%   CartesianState(4): vx                                           %
%   CartesianState(5): vy                                           %             
%   CartesianState(6): vz                                           %
%-------------------------------------------------------------------%
% Spherical states:                                                 % 
%   SphericalState(1): rho                                          %
%   SphericalState(2): theta                                        %
%   SphericalState(3): phi                                          %
%   SphericalState(4): d(rho)/dt = 0                                %
%   SphericalState(5): d(theta)/dt                                  %
%   SphericalState(6): d(phi)/dt                                    %
%-------------------------------------------------------------------%
if SphericalState6D(2) ~= 0
    error("No radial velocity. See Spherical2Cartesian.")
end
rho = SphericalState6D(1);
SphericalState4D = SphericalState6D([2,3,5,6]);

theta = SphericalState4D(1);
phi = SphericalState4D(2);

vTheta = SphericalState4D(3);
vPhi = SphericalState4D(4);
sTheta = sin(theta);
cTheta = cos(theta);
sPhi = sin(phi);
cPhi = cos(phi);

CartesianState = zeros(6, 1);

% Cartesian position and velocity
p = rho .* [cTheta .* sPhi;
            sTheta .* sPhi;
            cPhi];                          
v = rho .* [-vTheta .* sTheta .* sPhi + vPhi .* cTheta .* cPhi;         
            vTheta .* cTheta .* sPhi + vPhi .* sTheta .* cPhi;   
            -vPhi .* sPhi];
CartesianState(1:3) = p;
CartesianState(4:6) = v;

end