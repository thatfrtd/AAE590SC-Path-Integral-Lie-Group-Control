function L = cost(xtraj, xtarg, dt)
xf = xtraj(:, end);
cost_exit = norm(xf(1:2) - xtarg(1:2))^2; 

cost_t = sum(vecnorm(u_traj, 2, 1).^2) * dt;

% total cost
L = cost_exit + cost_t;
end
