function [f, B] = Euler_Poincare_matrices(G, J_b, B_map, u_ext)
%EULER_POINCARE Dynamics of twists (velocity in body frame)
%   The dynamical equation for the group in the body frame. Dual to the
%   Newton-Euler equations which are for dynamics in the space or inertial
%   frame.
%   G - group
%   twist - twist (velocity in body frame)
%   u - control (in body frame)
%   J_b - generalized inertia matrix
%   u_ext - known external control-like disturbances possibly sate dependent
arguments
    G
    J_b
    B_map
    u_ext = @(x) zeros([G.dim, 1])
end

% coadjoint action is always transpose of adjoint?
ad_star_twist = @(twist) G.ad(G, twist)';

%dtwist = J_b \ (ad_star_twist * J_b * twist + B * u);

f = @(x, twist) pagemtimes(J_b \ ad_star_twist(twist) * J_b, twist) + J_b \ u_ext(x);
B = @(t, x) J_b \ B_map;

end

