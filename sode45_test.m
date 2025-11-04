%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Travis Hastreiter 
% Created On: 3 November, 2025
% Description: Test of stochastic propagation without optimized control
% Most Recent Change: 3 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Brownian Motion Example
f = @(t, x, u) [0; 0] + u;
t_k = 0:0.01:1;
u_k = zeros([2, numel(t_k)]);

noise_delta_t = 1e-2;

sigma_1 = 2;
sigma_2 = 3;
sigma = [sigma_1, 0; 0, sigma_2];% * sqrt(delta_t);

x0 = [0; 0];

w = @(n) randn([2, n]);

m = 100;

x = zeros([2, numel(t_k), m]);
t = zeros([numel(t_k), m]);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

parfor i = 1:m
    [t(:, i), x(:, :, i)] = sode45(f, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances);
    i
end
%% Analyze distribution of trajectories to check if it matches expectations
figure
tiledlayout(1, 2);

nexttile
plot(t, squeeze(x(1, :, :))); hold on
plot(t, 3 * sigma_1 * sqrt(t), Color = "k"); hold on
plot(t, -3 * sigma_1 * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_1")
legend("x_1", "", "3 \sigma Bound")
title("x_1 Trajectories")
grid on

nexttile
plot(t, squeeze(x(2, :, :))); hold on
plot(t, 3 * sigma_2 * sqrt(t), Color = "k"); hold on
plot(t, -3 * sigma_2 * sqrt(t), Color = "k"); hold off
xlabel("Time")
ylabel("x_2")
legend("x_2", "", "3 \sigma Bound")
title("x_2 Trajectories")
grid on

%% Dynamical System Example
f = @(t, x, u) [x(2); u(1)];
t_k = 0:0.01:1;
u_k = -2 * ones([2, numel(t_k)]);

noise_delta_t = 1e-2;

w = @(n) randn([2, n]);

sigma_accel = 0.2; % [m / s2]
sigma = [0, sigma_accel; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

x0 = [0; 1];

m = 100;

x = zeros([2, numel(tspan), m]);
t = zeros([numel(tspan), m]);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

parfor i = 1:m
    [t(:, i), x(:, :, i)] = sode45(f, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances);
    i
end
%% Analyze distribution of trajectories to check if it matches expectations
figure
tiledlayout(1, 2);

nexttile
plot(t, squeeze(x(1, :, :))); hold on
xlabel("Time")
ylabel("Position")
legend("r")
title("Position vs Time")
grid on

nexttile
plot(t, squeeze(x(2, :, :))); hold on
xlabel("Time")
ylabel("Velocity")
legend("x_2")
title("Velocity vs Time")
grid on
