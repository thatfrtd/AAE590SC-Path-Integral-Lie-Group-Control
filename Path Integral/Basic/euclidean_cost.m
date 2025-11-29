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
    options.state_running_cost = @(x, v) 0
end

% Exit Cost
xf = xtraj(:, end, :);
vf = vtraj(:, end, :);
x_error = xf - xtarg;
v_error = vf - vtarg;
cost_exit = sqrt(squeeze(dot(x_error, pagemtimes(options.S_x, x_error)))) + sqrt(squeeze(dot(v_error, pagemtimes(options.S_v, v_error)))); 

% Control Running Cost
diff_cost = squeeze(sum(vecnorm(diff(u_traj, 1, 2) .^ 2), 2)); % penalize change in control 

cost_control_t = squeeze(sum(dot(u_traj, pagemtimes(R, u_traj))) * dt) + diff_cost * options.diff_multiplier;

% State Running Cost
cost_state_t = zeros(size(xtraj, 2:3));
for i = 1 : size(xtraj, 2)
    for j = 1 : size(xtraj, 3)
        cost_state_t(i, j) = options.state_running_cost(xtraj(:, i, j), vtraj(:, i, j));
    end
end
cost_state_t = sum(cost_state_t, 1)';

cost_t = cost_control_t + cost_state_t;

% total cost
L = 1 / 2 * cost_exit + 1 / 2 * cost_t;
end
