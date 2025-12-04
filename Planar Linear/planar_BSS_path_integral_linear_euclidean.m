%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 29 November, 2025
% Description: Test of euclidean belief sstate shield path integral control 
% on a planar double integrator system
% Most Recent Change: 29 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function u_k = planar_BSS_path_integral_linear_euclidean(integration_method, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, w, tolerances)

% integration_method = "oneeuler"; % "ode45" or "oneeuler"
% 
% %% Create Dynamics 
% z2 = zeros(2);
% I2 = eye(2);
% 
% A = [z2, I2; z2, z2];
% B = [z2; I2];
% 
% f = @(t, x, u) A * x + B * u;
% control_delta_t = 1e-2;
% t_k = 0:control_delta_t:1;
% u_k = 0 * ones([2, numel(t_k) - 1]);
% f_func = @(t, x) z2 * x;
% B_func = @(t, x) I2;
% 
% %% Define Initial Condition and Target
% x0 = [0; 0];
% v0 = [0; 1];
% xtarg = [1; 0]; 
% vtarg = [0; 0];

%% Define Noise
noise_delta_t = control_delta_t;%1e-2;

% w = @(n) randn([2, n]);

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


%% Time Histories
%sigma_accel = 0.5; % [m / s2]
%sigma = [sigma_accel, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

m = 50;

x = zeros([2, numel(t_k), m]);
v = zeros([2, numel(t_k), m]);
t = zeros([numel(t_k), m]);
delta_u = zeros([2, numel(t_k) - 1, m]);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);
for i = 1 : m
    x0_m_full = [x0; v0] + chol(P0, "lower") * randn([4, 1]); % Sample starting state

    [t(:, i), x_full(:, :), u_n(:, :, i), delta_u(:, :, i)] = sode45(f, u_k, sigma, w, t_k, noise_delta_t, x0_m_full, tolerances);
    x(:, :, i) = x_full(1:2, :);
    v(:, :, i) = x_full(3:4, :);
end
[t_nom, x2_nom, u_n, ~] = sode45(f, u_k, sigma*0, w, t_k, noise_delta_t, [x0; v0], tolerances);
x_nom = x2_nom(1:2, :);
v_nom = x2_nom(3:4, :);
%%
[t2, x2, v2, u_n2, delta_u2] = one_step_euler_maruyama_euclidean(f_func, B_func, u_k, sigma*0, t_k, x0, v0);

%%
figure
tiledlayout(1, 3);

nexttile
plot(t, squeeze(x(1, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
plot(t, squeeze(x(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off");
plot(t_nom, squeeze(x_nom(1, :)), Color = "r");
plot(t_nom, squeeze(x_nom(2, :)), Color = "b");
xlabel("Time")
ylabel("Position")
legend("r_x", "r_y")
title("Position vs Time")
grid on

nexttile
plot(t, squeeze(v(1, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
plot(t, squeeze(v(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off");
plot(t_nom, squeeze(v_nom(1, :)), Color = "r"); hold on
plot(t_nom, squeeze(v_nom(2, :)), Color = "b");
xlabel("Time")
ylabel("Velocity")
legend("v_x", "v_y")
title("Velocity vs Time")
grid on

nexttile
stairs(t_nom(1 : (end - 1)), squeeze(u_n(1, :) + delta_u(1, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
stairs(t_nom(1 : (end - 1)), squeeze(u_n(2, :) + delta_u(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off");
stairs(t_nom(1 : (end - 1)), squeeze(u_n(1, :)), Color = "r");
stairs(t_nom(1 : (end - 1)), squeeze(u_n(2, :)), Color = "b");
xlabel("Time")
ylabel("Control")
legend("a_x", "a_y")
title("Control Acceleration vs Time")
grid on

%% Trajectory
figure
tiledlayout(1, 2)

nexttile
plot(squeeze(x(1, :, :)), squeeze(x(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on

% Compute confidence ellipses
[P_eigvecs, P_eigvals] = pageeig(P_k(1:2, 1:2, :)); % Looks more correct then projecting ellipsoid...

% Create ellipses points
thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
ellipse_3sigma = zeros([2, 100, numel(t_k)]);
for k = 1:numel(t_k)
    ellipse_3sigma(:, :, k) = x_nom(1:2, k) + P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
    %X_k(:, :, k) = chol(P_k(:, :, k), "lower");
end

plot(squeeze(ellipse_3sigma(1, :, 2:end)), squeeze(ellipse_3sigma(2, :, 2:end)), Color = "k", HandleVisibility='off'); hold on
plot(squeeze(x_nom(1, :)), squeeze(x_nom(2, :)), Color = "r")
%plot(squeeze(x2(1, :)), squeeze(x2(2, :)), Color = "b", LineStyle="--")
plot(obstacle_radius * cos(0:0.01:2 * pi) + obstacle_center(1), obstacle_radius * sin(0:0.01:2 * pi) + obstacle_center(2), LineWidth=1,Color="k")



xlabel("X [m]")
ylabel("Y [m]")
title("Position Trajectory")
grid on
axis equal

nexttile
plot(squeeze(v(1, :, :)), squeeze(v(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
plot(squeeze(v_nom(1, :)), squeeze(v_nom(2, :)), Color = "r")
plot(squeeze(v2(1, :)), squeeze(v2(2, :)), Color = "b", LineStyle="--")
xlabel("X [m / s]")
ylabel("Y [m / s]")
title("Velocity Trajectory")
grid on
axis equal

%% Iterations
figure
plot(average_cost); hold on
scatter(1 : iterations, cost_t, MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
scatter(1 : iterations, cost_exit, MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
plot(mean(cost_t, 1), Color = "b");
plot(mean(cost_exit, 1), Color = "r");
grid on
legend("Total Cost", "Path Cost", "Exit Cost")
xlabel("Iteration")
ylabel("Cost")
title("Cost vs Iteration")
yscale("log")