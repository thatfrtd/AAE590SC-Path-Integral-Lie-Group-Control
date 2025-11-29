function [L, cost_t, cost_exit] = euclidean_cost(xtraj, vtraj, xtarg, vtarg, u_traj, dt, R, options)
arguments
    xtraj
    vtraj
    xtarg
    vtarg
    u_traj
    dt
    R
    options.S_x = 2 * eye(size(xtraj, 1))
    options.S_v = eye(size(vtraj, 1))
    options.diff_multiplier = 0
end
xf = xtraj(:, end, :);
vf = vtraj(:, end, :);
x_error = xf - xtarg;
v_error = vf - vtarg;
cost_exit = sqrt(squeeze(dot(x_error, pagemtimes(options.S_x, x_error)))) + sqrt(squeeze(dot(v_error, pagemtimes(options.S_v, v_error)))); 

diff_cost = squeeze(sum(vecnorm(diff(u_traj, 1, 2) .^ 2), 2)); % penalize change in control 

cost_t = squeeze(sum(dot(u_traj, pagemtimes(R, u_traj))) * dt) + diff_cost * options.diff_multiplier;

% total cost
L = 1 / 2 * cost_exit + 1 / 2 * cost_t;
end
