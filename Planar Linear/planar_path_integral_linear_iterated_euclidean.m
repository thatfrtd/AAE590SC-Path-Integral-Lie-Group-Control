%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 November, 2025
% Description: Test of euclidean path integral control on a planar double 
% integrator system
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [u_k, u_n, cost_t, cost_exit, average_cost] = planar_path_integral_linear_iterated_euclidean(integration_method, iterations, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, noise_delta_t, w, tolerances)

sigma_accel = 0.5; % [m / s2]
sigma = [sigma_accel, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

%% Run Monte Carlo
m = 50;

x = zeros([2, numel(t_k), m]);
v = zeros([2, numel(t_k), m]);
t = zeros([numel(t_k), m]);
u_n = zeros([2, numel(t_k) - 1]);
delta_u = zeros([2, numel(t_k) - 1, m]);

tolerances = odeset(RelTol=1e-4, AbsTol=1e-4, InitialStep=0.1, MaxStep=0.1);

%iterations = 210;
R = eye(2) * 1;
S_x = 800 * eye(2);
S_v = 200 * eye(2);
lambda_matrix = sigma * sigma' * R;
lambda = lambda_matrix(1)*100;
average_cost = zeros([1, iterations]);
cost_t = zeros([m, iterations]);
cost_exit = zeros([m, iterations]);
cost_state = zeros([numel(t_k) - 1, m]);
for j = 1 : iterations
    if integration_method == "ode45"
        for i = 1 : m
            [t(:, i), x_full(:, :), u_n(:, :, i), delta_u(:, :, i)] = sode45(f, u_k, sigma, w, t_k, noise_delta_t, [x0; v0], tolerances);
            x(:, :, i) = x_full(1:2, :);
            v(:, :, i) = x_full(3:4, :);
        end
    elseif integration_method == "oneeuler"
        parfor i = 1 : m
            [t(:, i), x(:, :, i), v(:, :, i), ~, delta_u(:, :, i)] = one_step_euler_maruyama_euclidean(f_func, B_func, u_k, sigma, t_k, x0, v0);
        end
    end
    
    [L, cost_t(:, j), cost_exit(:, j)] = iterative_euclidean_cost(x, v, xtarg, vtarg, u_k + delta_u, control_delta_t, R, S_x = S_x, S_v = S_v);

    %u_k = euclidean_update(u_k, delta_u, L, lambda);
    [u_k, L_new, L_raw] = iterative_euclidean_update(x, v, v * 0, B(3:4, :), u_k, delta_u, lambda, R, control_delta_t, cost_exit(:, j), cost_state);

    average_cost(j) = sum(L_raw,"all") / numel(L_raw);
    average_cost(j)

    if mod(j, 30) == 0 || j == 1
        [t_new, x_new, v_new, u_new, delta_u_new] = one_step_euler_maruyama_euclidean(f_func, B_func, u_k, sigma * 0, t_k, x0, v0);

        figure
        %plot(squeeze(x(1, :, :)), squeeze(x(2, :, :))); hold on
        colormap(jet);  
        %L_flat = repmat(L_new', 101, 1);
        L_new = [L_new; L_new(end, :)];
        scatter(squeeze(x(1, :)), squeeze(x(2, :)), [], L_new(:), "filled"); hold on
        alpha(L_new(:))
        plot(squeeze(x_new(1, :)), squeeze(x_new(2, :)), LineStyle="--", LineWidth=1, Color="r"); hold on
        scatter(xtarg(1), xtarg(2));
        xlabel("X [m]")
        ylabel("Y [m]")
        title("Trajectory")
        subtitle(sprintf("Iteration %d, avg cost %.3f", j, average_cost(j)))
        grid on
        axis equal
        colorbar
    end
end


end
