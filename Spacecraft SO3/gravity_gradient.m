function [moment] = gravity_gradient(x, orbit_mean_motion, I_b)
%GRAVITY_GRADIENT Create gravity gradient moment for circular orbit
%   Calculate the moment from the gravity gradient on a spacecraft in a
%   circular orbit. Assumes the rotation matrix is the orientation of the
%   spacecraft with respect to the orbit frame so body frame to orbit
%   frame.
arguments (Input)
    x % Any element of a Lie group that has a 3D rotation matrix (DCM)
    orbit_mean_motion
    I_b
end

% Get DCM
R = x.R'; % Transpose because equations expect orbit to body DCM

% Orbit radial vector
ohat_1 = R(:, 1);

% Create gravity gradient moment
moment = 3 * orbit_mean_motion ^ 2 * cross(ohat_1, I_b * ohat_1);
end