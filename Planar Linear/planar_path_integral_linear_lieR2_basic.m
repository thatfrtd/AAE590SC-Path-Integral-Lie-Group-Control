%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 November, 2025
% Description: Test of lie group path integral control on a planar double 
% integrator system represented by the Lie Group R2
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function u_k = planar_path_integral_linear_lieR2_basic(z2, I2, tolerances, noise_delta_t, w, f, B, f_with_control, control_delta_t, t_k, u_k, x0, xtarg, g0, twist0, gtarg, twisttarg)

% z2 = zeros(2);
% I2 = eye(2);
% 
% tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);
% 
% %% Define Noise
% noise_delta_t = 1e-2;
% 
% w = @(n) randn([2, n]);

sigma_accel = 0.5; % [kg m / s2]
sigma = [sigma_accel, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

% %% Create Dynamics 
% mass = 1;
% J_b = mass; % Generalized inertia
% 
% B_map = I2; 
% [f, B] = Euler_Poincare_matrices(R2(), J_b, B_map);
% 
% f_with_control = @(t, x, u) [z2, I2; z2, z2] * x + [z2; B(t, x)] * u;
% control_delta_t = 0.01;
% t_k = 0:control_delta_t:1;
% u_k = 0 * ones([2, numel(t_k) - 1]);
% 
% %% Define Initial Condition and Target
% x0 = [0; 0; 0; 1];
% xtarg = [1; 0; 0; 0];
% g0 = R2(x0(1:2));
% twist0 = x0(3:4);

% gtarg = R2(xtarg(1:2));
% twisttarg = xtarg(3:4);

%% Run Monte Carlo
m = 50;

g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
twist_k = zeros(g0.dim, numel(t_k), m); % array of twists (body velocities)
delta_u = zeros([2, numel(t_k) - 1, m]);
w_k = zeros([2, numel(t_k) - 1, m]);

iterations = 250;
R = eye(2);
S_g = 1000 * eye(2);
S_twist = 250 * eye(2);
lambda_matrix = sigma * sigma' * R;
lambda = lambda_matrix(1) * 100;
average_cost = zeros([1, iterations]);
cost_exit = zeros([m, iterations]);
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

%% Time Histories
sigma_accel = 0.2; % [m / s2]
sigma = [sigma_accel, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

m = 50;

x = zeros([4, numel(t_k), m]);
t = zeros([numel(t_k), m]);
delta_u = zeros([2, numel(t_k) - 1, m]);

for i = 1 : m
    [t(:, i), x(:, :, i), ~, delta_u(:, :, i)] = sode45(f_with_control, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances);
end
[t_nom, x_nom, u_n, ~] = sode45(f_with_control, u_k, sigma*0, w, t_k, noise_delta_t, x0, tolerances);

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
plot(t, squeeze(x(3, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
plot(t, squeeze(x(4, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off");
plot(t_nom, squeeze(x_nom(3, :)), Color = "r"); hold on
plot(t_nom, squeeze(x_nom(4, :)), Color = "b");
xlabel("Time")
ylabel("Velocity")
legend("v_x", "v_y")
title("Velocity vs Time")
grid on

nexttile
stairs(t(1 : end - 1, :), squeeze(u_n(1, :) + delta_u(1, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
stairs(t(1 : end - 1, :), squeeze(u_n(2, :) + delta_u(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off");
stairs(t_nom(1 : end - 1, :), squeeze(u_n(1, :)), Color = "r");
stairs(t_nom(1 : end - 1, :), squeeze(u_n(2, :)), Color = "b");
xlabel("Time")
ylabel("Control")
legend("a_x", "a_y")
title("Control Acceleration vs Time")
grid on

%% Trajectory
[g_nom, twist_k_nom] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma * 0);
x_g_nom = [g_nom.element];

figure
plot(squeeze(x(1, :, :)), squeeze(x(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
plot(squeeze(x_nom(1, :)), squeeze(x_nom(2, :)), Color = "r")
plot(squeeze(x_g_nom(1, :)), squeeze(x_g_nom(2, :)), Color = "g")
xlabel("X [m]")
ylabel("Y [m]")
title("Trajectory")
grid on
axis equal

%% Iterations
figure
plot(average_cost); hold on
scatter(1 : iterations, squeeze(sum(cost_t, 1)), MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
scatter(1 : iterations, cost_exit, MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
plot(mean(squeeze(sum(cost_t, 1)), 1), Color = "b");
plot(mean(cost_exit, 1), Color = "r");
grid on
legend("Total Cost", "Path Cost", "Exit Cost")
xlabel("Iteration")
ylabel("Cost")
title("Cost vs Iteration")
yscale("log")
end