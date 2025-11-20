function [f, B] = Euler_Poincare_matrices(G, J_b, B_map)
%EULER_POINCARE Dynamics of twists (velocity in body frame)
%   The dynamical equation for the group in the body frame. Dual to the
%   Newton-Euler equations which are for dynamics in the space or inertial
%   frame.
%   G - group
%   twist - twist (velocity in body frame)
%   u - control (in body frame)
%   J_b - generalized inertia matrix


% coadjoint action is always transpose of adjoint?
ad_star_twist = @(twist) G.ad(twist)';

%dtwist = J_b \ (ad_star_twist * J_b * twist + B * u);

f = @(x, twist) pagemtimes(J_b \ ad_star_twist(twist) * J_b, twist);
B = @(x, twist) J_b \ B_map;

end

