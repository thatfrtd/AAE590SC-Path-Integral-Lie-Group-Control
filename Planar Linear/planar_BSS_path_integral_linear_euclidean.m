%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 29 November, 2025
% Description: Test of euclidean belief sstate shield path integral control 
% on a planar double integrator system
% Most Recent Change: 29 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u_k, u_n, cost_t, cost_exit, average_cost] = planar_BSS_path_integral_linear_euclidean(integration_method, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, w, tolerances)

%% Define Noise
noise_delta_t = control_delta_t;%1e-2;

sigma_accel = 0.3; % [m / s2]
sigma = [sigma_accel*1.5, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

%% Create Kalman Filter - completely deterministic because of linear system
% Define initial state covariance
sigma_xhat0 = [10e-3; ... % r_x
            2e-3; ... % r_y
            6e-3; ... % v_x
            3e-3]; ... % v_y
P0 = diag(sigma_xhat0 .^ 2);

% Discretize system - going to do rough approximation (should be actual discretization)

P_k = zeros([4, 4, numel(t_k)]);
P_k(:, :, 1) = P0;

filter.A_k = eye(4) + A * control_delta_t;
filter.G_k = B() * sigma / sqrt(noise_delta_t) * control_delta_t;
for k = 1 : numel(t_k) - 1
    P_k(:, :, k + 1) = filter.A_k * P_k(:, :, k) * filter.A_k' ...
                          + filter.G_k * filter.G_k';
end

% filter.B_k = B() * control_delta_t;
% filter.c_k = zeros([4, 1]);
% filter.C_k = eye([4, 4]); % Full state measurement
% filter.D_k = diag([5e-1, 5e-1, 5e-2, 5e-2]);

% Compute Kalman filter covariance
% [Ptilde_minus_k, P_k, filter.L_k] = compute_Kalman_matrices_apriori(P0, filter.A_k, filter.G_k, filter.C_k, filter.D_k, numel(t_k));

%% Define Obstacles
obstacle_center = [0.6; 0.1];
obstacle_radius = 0.3;
obstacle_penalty = 1;
%logistic = @(x) 1e4 ./ (1 + exp(-1000*x)); % Smooth out penalty
logistic = @(x) 1 ./ (1 + exp(-8000*x)); % Smooth out penalty
obstacle_constraint = @(x) obstacle_radius - sqrt(dot(x - obstacle_center, x - obstacle_center)); % < 0
state_running_cost = @(x, v) obstacle_penalty * logistic(obstacle_radius - sqrt(dot(x - obstacle_center, x - obstacle_center)));
C = 1;

% Calculate safe set gradient
state_sym = sym("state", [4, 1]);
state_running_cost_grad = matlabFunction(jacobian(state_running_cost(state_sym(1:2), state_sym(3:4)), state_sym), "Vars", {state_sym});
state_running_cost_hess = matlabFunction(hessian(state_running_cost(state_sym(1:2), state_sym(3:4)), state_sym), "Vars", {state_sym});

% Define backoff coefficient
P_fail = 0.01;
backoff_coef = sqrt(chi2inv(1 - P_fail, 2));

beta = 0.5;

% Define discrete control barrier function heuristic
h_func = chance_constraint_heuristic(@(x) state_running_cost(x(1:2), x(3:4)), @(x) state_running_cost_grad(x)', state_running_cost_hess, backoff_coef);
chance_constrained_state_running_cost = @(x, v, P) obstacle_penalty * logistic(obstacle_constraint(x) + backoff_coef * sqrt((x - obstacle_center)' * P(1:2, 1:2) * (x - obstacle_center) / norm(x - obstacle_center) ^ 2));

%% Run Monte Carlo
m = 100;

x = zeros([2, numel(t_k), m]);
v = zeros([2, numel(t_k), m]);
t = zeros([numel(t_k), m]);
u_n = zeros([2, numel(t_k) - 1, m]);
delta_u = zeros([2, numel(t_k) - 1, m]);

tolerances = odeset(RelTol=1e-4, AbsTol=1e-4, InitialStep=0.1, MaxStep=0.1);

iterations = 400;
R = eye(2) * 0.5;
S_x = 6000 * eye(2);
S_v = 1500 * eye(2);
lambda_matrix = sigma * sigma' * R;
lambda = lambda_matrix(1)*100;
average_cost = zeros([1, iterations]);
cost_t = zeros([m, iterations]);
cost_exit = zeros([m, iterations]);

% Copy covariance matrix for each rollout
P_k = repmat(P_k, 1, 1, 1, m);

%%
for j = 1 : iterations
    if integration_method == "ode45"
        for i = 1 : m
            [t(:, i), x_full(:, :), u_n(:, :, i), delta_u(:, :, i)] = sode45(f, u_k, sigma, w, t_k, noise_delta_t, [x0; v0], tolerances);
            x(:, :, i) = x_full(1:2, :);
            v(:, :, i) = x_full(3:4, :);
        end
    elseif integration_method == "oneeuler"
        parfor i = 1 : m
            %x0_m_full = [x0; v0] + chol(P0, "lower") * randn([4, 1]); % Sample starting state

            [t(:, i), x(:, :, i), v(:, :, i), u_n(:, :, i), delta_u(:, :, i)] = one_step_euler_maruyama_euclidean(f_func, B_func, u_k, sigma, t_k, x0, v0);

            % Kalman filter is deterministic so same for each rollout
            %[xtilde_k, P_k(:, :, :, i)] = get_state_update_Kalman(x(:, :, i), P0, u_k, filter);
        end
    end
   
    [L, cost_t(:, j), cost_exit(:, j)] = euclidean_cost(x, v, xtarg, vtarg, u_n + delta_u, control_delta_t, R, S_x = S_x, S_v = S_v, state_running_cost = chance_constrained_state_running_cost, P_k = P_k);
    %[L, cost_t(:, j), cost_exit(:, j)] = euclidean_cost(x, v, xtarg, vtarg, u_n + delta_u, control_delta_t, R, S_x = S_x, S_v = S_v, state_running_cost = state_running_cost);
    
    %C_safe = safe_condition_violation_cost(C, h_func, beta, [x; v], P_k);
    %L = L + sum(C_safe, 1)'; % Add discrete control barrier function heuristic cost

    u_k = euclidean_update(u_k, delta_u, L, lambda);

    average_cost(j) = sum(L) / numel(L);
    average_cost(j)

    if mod(j, 30) == 0 || j == 1
        [t_new, x_new, v_new, u_new, delta_u_new] = one_step_euler_maruyama_euclidean(f_func, B_func, u_k, sigma * 0, t_k, x0, v0);

        figure
        %plot(squeeze(x(1, :, :)), squeeze(x(2, :, :))); hold on
        colormap(jet);  
        L_flat = repmat(L', numel(t_k), 1);
        scatter(squeeze(x(1, :)), squeeze(x(2, :)), [], L_flat(:), "filled"); hold on
        plot(squeeze(x_new(1, :)), squeeze(x_new(2, :)), LineStyle="--", LineWidth=1, Color="r"); hold on
        plot(obstacle_radius * cos(0:0.01:2 * pi) + obstacle_center(1), obstacle_radius * sin(0:0.01:2 * pi) + obstacle_center(2), LineWidth=1,Color="k")
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