%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 November, 2025
% Description: Test of lie group path integral control on a planar double 
% integrator system represented by the Lie Group R2
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% function u_k = planar_path_integral_linear_lieR2(z2, I2, tolerances, noise_delta_t, w, f, B, f_with_control, control_delta_t, t_k, u_k, x0, xtarg, g0, twist0, gtarg, twisttarg)

function [u_k, cost_t, cost_exit, average_cost]  = planar_path_integral_linear_lieR2(z2, I2, tolerances, noise_delta_t, w, f, B, f_with_control, control_delta_t, t_k, u_k, xtarg, g0, twist0, gtarg, twisttarg, iterations)

sigma_accel = 0.5; % [kg m / s2]
sigma = [sigma_accel, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

x0 = [g0.element; twist0];
%% Run Monte Carlo
m = 50;

g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
twist_k = zeros(g0.dim, numel(t_k), m); % array of twists (body velocities)
delta_u = zeros([2, numel(t_k) - 1, m]);
w_k = zeros([2, numel(t_k) - 1, m]);

%iterations = 250;
R = eye(2);
S_g = 8000 * eye(2);
S_twist = 2000 * eye(2);
lambda_matrix = sigma * sigma' * R;
lambda = lambda_matrix(1)*100;
average_cost = zeros([1, iterations]);
cost_exit = zeros([m, iterations]);
%%Different Cost_t 
cost_t = zeros([numel(t_k) - 1, m, iterations]);

for j = 1 : iterations
    parfor i = 1:m
        % [t_k2, x_k2, v_k2, u_k2, delta_u2, w_k] = one_step_euler_maruyama_euclidean(@(t,x)0, B, u_k, sigma, t_k, x0(1:2), x0(3:4));

        [g_k(:, i), twist_k(:, :, i), delta_u(:, :, i), w_k(:, :, i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
    end
    
    % for i = 1:m
    %     [g_k(:, i), twist_k(:, :, i), delta_u(:, :, i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
    % 
    % end

    [L, cost_t(:, :, j), cost_exit(:, j)] = liegroup_cost(g_k, twist_k, u_k + delta_u, control_delta_t, gtarg, twisttarg, R, S_g = S_g, S_twist = S_twist);

%    [u_k, D_k, L_k] = liegroup_update(g_k, u_k, twist_k, pagemtimes(sigma, w_k), f, B, R, L, t_k, lambda);
    [u_k, D_k, L_k] = liegroup_update(g_k, u_k, twist_k, delta_u, f, B, R, L, t_k, lambda);

    average_cost(j) = sum(L_k, "all") / numel(L_k);
    average_cost(j)

    
end

[t_opt, x_opt, ~, delta_u_opt] = sode45(f_with_control, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances);

      

figure
x = [g_k(2:end, :).element];
colormap(jet);  
scatter(squeeze(x(1, :)), squeeze(x(2, :)), [], D_k(:), "filled"); hold on
alpha(D_k(:));
plot(x_opt(1, :), x_opt(2, :));
%scatter(xtarg(1), xtarg(2));
xlabel("X [m]")
ylabel("Y [m]")
title("Trajectory")
subtitle(sprintf("Iteration %d, avg cost %.3f", j, average_cost(j)))
grid on
axis equal
colorbar
end
