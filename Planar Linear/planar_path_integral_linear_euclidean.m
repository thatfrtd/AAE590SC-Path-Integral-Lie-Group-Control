%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 November, 2025
% Description: Test of euclidean path integral control on a planar double 
% integrator system
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

integration_method = "oneeuler"; % "ode45" or "oneeuler"

%% Create Dynamics 
z2 = zeros(2);
I2 = eye(2);

A = [z2, I2; z2, z2];
B = [z2; I2];

f = @(t, x, u) A * x + B * u;
control_delta_t = 0.01;
t_k = 0:control_delta_t:1;
u_k = 0 * ones([2, numel(t_k) - 1]);
f_func = @(t, x) z2 * x;
B_func = @(t, x) I2;

%% Define Initial Condition and Target
x0 = [0; 0];
v0 = [0; 0.1];
xtarg = [1; 0]; 
vtarg = [0; 0];

%% Define Noise
noise_delta_t = 1e-2;

w = @(n) randn([2, n]);

sigma_accel = 0.4; % [m / s2]
sigma = [sigma_accel, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

%% Run Monte Carlo
m = 200;

x = zeros([2, numel(t_k), m]);
v = zeros([2, numel(t_k), m]);
t = zeros([numel(t_k), m]);
u_n = zeros([2, numel(t_k) - 1, m]);
delta_u = zeros([2, numel(t_k) - 1, m]);

tolerances = odeset(RelTol=1e-4, AbsTol=1e-4, InitialStep=0.1, MaxStep=0.1);

iterations = 200;
lambda = 0.01;
eta = 0.9; % 0 means only path cost (don't do), 1 means only terminal cost
average_cost = zeros([1, iterations]);
cost_t = zeros([m, iterations]);
cost_exit = zeros([m, iterations]);

for j = 1 : iterations
    if integration_method == "ode45"
        for i = 1 : m
            [t(:, i), x_full(:, :), u_n(:, :, i), delta_u(:, :, i)] = sode45(f, u_k, sigma, w, t_k, noise_delta_t, [x0; v0], tolerances);
            x(:, :, i) = x_full(1:2, :);
            v(:, :, i) = x_full(3:4, :);
        end
    elseif integration_method == "oneeuler"
        parfor i = 1 : m
            [t(:, i), x(:, :, i), v(:, :, i), u_n(:, :, i), delta_u(:, :, i)] = one_step_euler_maruyama_euclidean(f_func, B_func, u_k, sigma, t_k, x0, v0);
        end
    end
    
    [L, cost_t(:, j), cost_exit(:, j)] = euclidean_cost(x, v, xtarg, vtarg, u_n + delta_u, control_delta_t, eta);
    u_k = euclidean_update(u_k, delta_u, L, lambda);

    average_cost(j) = sum(L) / numel(L);
    average_cost(j)

    if mod(j, 50) == 0 || j == 1
        figure
        plot(squeeze(x(1, :, :)), squeeze(x(2, :, :))); hold on
        scatter(xtarg(1), xtarg(2));
        xlabel("X [m]")
        ylabel("Y [m]")
        title("Trajectory")
        subtitle(sprintf("Iteration %d, avg cost %.3f", j, average_cost(j)))
        grid on
        axis equal
    end
end


%% Time Histories
sigma_accel = 0.2; % [m / s2]
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
[t2, x2, v2, u_n2, delta_u2] = one_step_euler_maruyama_euclidean(f_func, B_func, u_k, sigma, t_k, x0, v0);

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
