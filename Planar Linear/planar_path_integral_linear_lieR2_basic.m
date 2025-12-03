%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 November, 2025
% Description: Test of lie group path integral control on a planar double 
% integrator system represented by the Lie Group R2
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u_k, average_cost, cost_t, cost_exit, g_k, twist_k, D_k, L_k] = ...
    planar_path_integral_linear_lieR2_basic( ...
    f, B, g0, twist0, u_k, t_k, sigma, gtarg, twisttarg, control_delta_t, ...
    R, S_g, S_twist, lambda, iterations, m, f_with_control, w, noise_delta_t, ...
        x0, tolerances)

%%Different Cont_t
%cost_t = zeros([numel(t_k) - 1, m, iterations]);
cost_t = zeros([m, 1, iterations]);

for j = 1 : iterations
    parfor i = 1 : m
        % [t_k2, x_k2, v_k2, u_k2, delta_u2, w_k] = one_step_euler_maruyama_euclidean(@(t,x)0, B, u_k, sigma, t_k, x0(1:2), x0(3:4));

        [g_k(:, i), twist_k(:, :, i), delta_u(:, :, i), w_k(:,:,i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
    end
    
    %[L, cost_t(:, :, j), cost_exit(:, j)] = liegroup_cost(g_k, twist_k, u_k + delta_u, control_delta_t, gtarg, twisttarg, R, S_g = S_g, S_twist = S_twist);
    [L_k, cost_t(:, :, j), cost_exit(:, j)] = basic_Lie_cost(g_k, twist_k, u_k + delta_u, control_delta_t, gtarg, twisttarg, R, S_g = S_g, S_twist = S_twist);

    %[u_k, D_k, L_k] = liegroup_update(g_k, u_k, twist_k, pagemtimes(sigma, w_k), f, B, R, L, t_k, lambda);
    [u_k, D_k] = euclidean_update(u_k, delta_u, L_k, lambda);
    D_k = repmat(D_k', numel(t_k) - 1, 1);

    average_cost(j) = sum(L_k, "all") / numel(L_k);
    average_cost(j)

    if mod(j, 30) == 0 || j == 1
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
end

