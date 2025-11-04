function [L, cost_t, cost_exit] = euclidean_cost(xtraj, xtarg, u_traj, dt, eta, diff_multiplier)
arguments
    xtraj
    xtarg
    u_traj
    dt
    eta = 0.9
    diff_multiplier = 0.02;
end
xf = xtraj(:, end, :);
cost_exit = squeeze(vecnorm(xf - xtarg)); 

diff_cost = squeeze(sum(vecnorm(diff(u_traj, 1, 2)), 2)); % penalize change in control 

cost_t = squeeze(sum(vecnorm(u_traj, 2, 1)) * dt) + diff_cost * diff_multiplier;

% total cost
L = cost_exit * eta + cost_t * (1 - eta);
end
