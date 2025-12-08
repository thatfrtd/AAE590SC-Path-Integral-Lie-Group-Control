%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 December, 2025
% Description: Main script to run the Planar Path Integral Functions
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z2 = zeros(2);
I2 = eye(2);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

noise_delta_t = 1e-2;

w = @(n) randn([2, n]);

mass = 1;
J_b = mass;

B_map = I2; 
[f, B] = Euler_Poincare_matrices(R2(), J_b, B_map);

f_with_control = @(t, x, u) [z2, I2; z2, z2] * x + [z2; B(t, x)] * u;
control_delta_t = 0.01;
t_k = 0:control_delta_t:1;
%disp(size(t_k))
u_k = 0 * ones([2, numel(t_k) - 1]);

x0 = [0; 0; 0; 1];
xtarg = [1; 0; 0; 0];
g0 = R2(x0(1:2));
twist0 = x0(3:4);

gtarg = R2(xtarg(1:2));
twisttarg = xtarg(3:4);

iterations = 2;

max_steps = 50;

for i = 1 : max_steps
    % Update Control
    [u_k, cost_t, cost_exit, average_cost]  = planar_path_integral_linear_lieR2(z2, I2, tolerances, noise_delta_t, w, f, B, f_with_control, control_delta_t, t_k, u_k, xtarg, g0, twist0, gtarg, twisttarg, iterations);


    % Step forward one timestep
    [g_k_next, twist_k_next, ~, ~] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, [0, control_delta_t], sigma);
    g0 = g_k_next(:, 2);
    twist0 = twist_k_next(:, 2); 
    u_k = [u_k(:,2:end), zeros(size(u_k,1), 1)];

    % Check if stopping condition met
    g_error = gtarg.left_invariant_error(g0);
    twist_error = twisttarg - twist0;
end

%[u_k, cost_t, cost_exit, average_cost]  = planar_path_integral_linear_lieR2_basic(z2, I2, tolerances, noise_delta_t, w, f, B, f_with_control, control_delta_t, t_k, u_k, x0, xtarg, g0, twist0, gtarg, twisttarg, iterations);

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