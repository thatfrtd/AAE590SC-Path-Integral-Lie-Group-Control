function [L, cost_t, cost_exit] = liegroup_cost(g_k, twist_k, u_k, delta_t, g_desired, twist_desired, eta, diff_multiplier)
%LIEGROUP_COST Summary of this function goes here
%   Cost of a trajectory on a lie group, should be valid for any lie group
arguments
    g_k
    twist_k
    u_k
    delta_t
    g_desired
    twist_desired
    eta = 0.9
    diff_multiplier = 0.02;
end

%% Path Cost
diff_cost = squeeze(sum(vecnorm(diff(u_k, 1, 2)), 2)); % penalize change in control 

cost_t = squeeze(sum(vecnorm(u_k, 2, 1)) * delta_t) + diff_cost * diff_multiplier;

%% Exit Cost
g_f = g_k(end, :);
group_error = vecnorm(g_desired.left_invariant_error(g_f)); % when use right invariant?
twist_error = vecnorm(twist_desired - twist_k);

cost_exit = group_error + twist_error;

%% Total Cost
L = cost_exit * eta + cost_t * (1 - eta);

end