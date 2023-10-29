%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Cartesian coordinates ---> Spherical coordinates
% Cartesian states:
%   CartesianState(1): x
%   CartesianState(2): y
%   CartesianState(3): z
%   CartesianState(4): vx
%   CartesianState(5): vy
%   CartesianState(6): vz
% Spherical states:
%   SphericalState(1): rho
%   SphericalState(2): theta
%   SphericalState(3): phi
%   SphericalState(4): d(rho)/dt
%   SphericalState(5): d(theta)/dt
%   SphericalState(6): d(phi)/dt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SphericalState6D] = Cartesian2Spherical(CartesianState)
SphericalState6D = zeros(6,1);

SphericalState6D(1) = sqrt(CartesianState(1)^2 ...
                        + CartesianState(2)^2 ...
                        + CartesianState(3)^2);
SphericalState6D(2) = atan2(CartesianState(2), CartesianState(1));

SphericalState6D(3) = acos(CartesianState(3) / SphericalState6D(1));

SphericalState6D(5) = -sin(SphericalState6D(2)) / sin(SphericalState6D(3)) * CartesianState(4) ...
                    + cos(SphericalState6D(2)) / sin(SphericalState6D(3)) * CartesianState(5);

SphericalState6D(6) = cos(SphericalState6D(2)) * cos(SphericalState6D(3)) * CartesianState(4) ...
                    + sin(SphericalState6D(2)) * cos(SphericalState6D(3)) * CartesianState(5) ...
                    - sin(CartesianState(3)) * CartesianState(6);
end