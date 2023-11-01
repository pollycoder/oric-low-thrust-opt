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
%   SphericalState(4): vRho
%   SphericalState(5): vTheta
%   SphericalState(6): vPhi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [SphericalState6D] = Cartesian2Spherical(CartesianState, tol)
if nargin < 2
    tol = 1e-6;
end

SphericalState6D = zeros(6,1);

x1 = CartesianState(1);
x2 = CartesianState(2);
x3 = CartesianState(3);
x4 = CartesianState(4);
x5 = CartesianState(5);
x6 = CartesianState(6);

% rho
SphericalState6D(1) = sqrt(x1^2 + x2^2 + x3^2);
rho = SphericalState6D(1);

% theta
SphericalState6D(2) = atan2(x2, x1);
theta = mod(SphericalState6D(2), 2 * pi);

% phi
SphericalState6D(3) = acos(x3 / rho);
phi = mod(SphericalState6D(3), 2 * pi);

if abs(sin(phi)) < tol
    error("The Phi you chose makes the problem singular.");
end

% Trigonometric function
sTheta = sin(theta);
cTheta = cos(theta);
sPhi = sin(phi);
cPhi = cos(phi);


% vRho = 0
% vTheta
SphericalState6D(5) = -sTheta / sPhi * x4 + cTheta / sPhi * x5;

% vPhi
SphericalState6D(6) = cTheta * cPhi * x4 + sTheta * cPhi * x5 - sPhi * x6;

end