function L = cost(xtraj, xtarg, u_traj, dt)
xf = xtraj(:, end, :);
cost_exit = squeeze(vecnorm(xf - xtarg).^2); 

cost_t = squeeze(sum(vecnorm(u_traj, 2, 1).^2) * dt);

% total cost
L = cost_exit + cost_t * 0;
end
