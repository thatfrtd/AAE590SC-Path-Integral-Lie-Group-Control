function [L, cost_t, cost_exit] = liegroup_cost(g_k, twist_k, u_k, delta_t, g_desired, twist_desired, R, options)
%LIEGROUP_COST Summary of this function goes here
%   Cost of a trajectory on a lie group, should be valid for any lie group
arguments
    g_k
    twist_k
    u_k
    delta_t
    g_desired
    twist_desired
    R
    options.S_g = 500 * eye(size(twist_k, 1))
    options.S_twist = 100 * eye(size(twist_k, 1))
end

%% Path Cost
%diff_cost = squeeze(cumsum(vecnorm(diff(u_k, 1, 2)), 2)); % penalize change in control 

cost_t = 1 / 2 * cumsum(squeeze(dot(u_k, pagemtimes(R, u_k)) * delta_t), 1, "reverse");% + diff_cost * diff_multiplier;

%% Exit Cost
g_f = g_k(end, :);
group_error_cost = zeros([numel(g_f), 1]);
for i = 1 : numel(g_f)
    g_error = g_desired.left_invariant_error(g_f(i));
    group_error_cost(i) = 1 / 2 * sqrt(dot(g_error, options.S_g * g_error)); % when use right invariant?
end
twist_error = twist_desired - twist_k(:, end, :);
twist_error_cost = 1 / 2 * sqrt(squeeze(dot(twist_error, pagemtimes(options.S_twist, twist_error))));

cost_exit = (group_error_cost + twist_error_cost)';

%% Total Cost
L = cost_exit + cost_t;

end