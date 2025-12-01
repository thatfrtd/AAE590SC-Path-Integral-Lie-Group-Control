function [L, cost_t, cost_exit] = basic_Lie_cost(g_k, twist_k, u_traj, dt, g_desired, twist_desired, R, options)
arguments
    g_k
    twist_k
    u_traj
    dt
    g_desired
    twist_desired
    R
    options.S_g = 2 * eye(g_k(1).dim)
    options.S_twist = eye(size(twist_k, 1))
    options.diff_multiplier = 0
end

g_f = g_k(end, :);
group_error_cost = zeros([numel(g_f), 1]);
for i = 1 : numel(g_f)
    g_error = g_desired.left_invariant_error(g_f(i));
    group_error_cost(i) = 1 / 2 * sqrt(dot(g_error, options.S_g * g_error)); % when use right invariant?
end
twist_error = twist_desired - twist_k(:, end, :);
twist_error_cost = 1 / 2 * sqrt(squeeze(dot(twist_error, pagemtimes(options.S_twist, twist_error))));

cost_exit = (group_error_cost + twist_error_cost);

diff_cost = squeeze(sum(vecnorm(diff(u_traj, 1, 2) .^ 2), 2)); % penalize change in control 

cost_t = 1 / 2 * squeeze(sum(dot(u_traj, pagemtimes(R, u_traj))) * dt) + diff_cost * options.diff_multiplier;

% total cost
L = cost_exit + cost_t;
end
