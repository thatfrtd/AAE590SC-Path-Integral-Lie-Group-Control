%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 December, 2025
% Description: Main script to iterating and planar PI eucliden 
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%integration_method = "oneeuler";

z2 = zeros(2);
I2 = eye(2);

integration_method = "oneeuler"; % "ode45" or "oneeuler"

A = [z2, I2; z2, z2];
B = [z2; I2];

f = @(t, x, u) pagemtimes(A, x) + pagemtimes(B, u);
control_delta_t = 0.01;
t_k = 0:control_delta_t:1;
u_k = 0 * ones([2, numel(t_k) - 1]);
f_func = @(t, x) z2 * x;
B_func = @(t, x) I2;

x0 = [0; 0];
v0 = [0; 1];
xtarg = [1; 0]; 
vtarg = [0; 0];

noise_delta_t = 1e-2;

w = @(n) randn([2, n]);

tolerances = odeset(RelTol=1e-4, AbsTol=1e-4, InitialStep=0.1, MaxStep=0.1);

iterations = 400;

%[u_k, u_n, cost_t, cost_exit, average_cost]  = planar_path_integral_linear_iterated_euclidean(integration_method, iterations, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, noise_delta_t, w, tolerances);

%[u_k, u_n, cost_t, cost_exit, average_cost] = planar_path_integral_linear_euclidean(integration_method, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, noise_delta_t, w, tolerances);
% 
[u_k, u_n, cost_t, cost_exit, average_cost] = planar_BSS_path_integral_linear_euclidean(integration_method, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, w, tolerances);
% 
% %% Time Histories
sigma_accel = 0.2; %0.3 for PI euclidean 
sigma = [sigma_accel, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

m = 50;

x = zeros([2, numel(t_k), m]);
v = zeros([2, numel(t_k), m]);
t = zeros([numel(t_k), m]);
delta_u = zeros([2, numel(t_k) - 1, m]);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

for i = 1 : m
    [t(:, i), x_full(:, :), u_n(:, :, i), delta_u(:, :, i)] = sode45(f, u_k, sigma, w, t_k, noise_delta_t, [x0; v0], tolerances);
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
plot(squeeze(x_nom(1, :)), squeeze(x_nom(2, :)), Color = "r")
plot(squeeze(x2(1, :)), squeeze(x2(2, :)), Color = "b", LineStyle="--")
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
hold on
scatter(1 : iterations, cost_t, MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
scatter(1 : iterations, cost_exit, MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
plot(mean(cost_t, 1), Color = "b");
plot(mean(cost_exit, 1), Color = "r");
plot(average_cost, Color = "m"); 
grid on
legend("Path Cost", "Exit Cost", "Average Cost")
xlabel("Iteration")
ylabel("Cost")
title("Cost vs Iteration")
yscale("log")

