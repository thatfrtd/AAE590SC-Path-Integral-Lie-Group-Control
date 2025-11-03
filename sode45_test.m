%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590ACA
% Stochastic SCP Rocket Landing Project
% Author: Travis Hastreiter 
% Created On: 15 April, 2025
% Description: Test of stochastic propagation without optimized feedback 
% control
% Most Recent Change: 15 April, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Brownian Motion Example
f = @(t, x, u, p) [0; 0];
u = @(t, x) [0; 0];
p = 0;

delta_t = 1e-3;

sigma_1 = 2;
sigma_2 = 3;
G = @(t, x, u, p) [sigma_1, 0; 0, sigma_2];% * sqrt(delta_t);

x0 = [0; 0];

w = @(n) randn([2, n]);

tspan = 0:delta_t:1;

m = 100;

x = zeros([2, numel(tspan), m]);
t = zeros([numel(tspan), m]);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

t_k = tspan;

parfor i = 1:m
    w_k = w(numel(t_k));
    w_func = @(t) w_k(:, floor(t /delta_t) + 1);
    [t(:, i), x(:, :, i)] = sode45(f, G, u, p, w, tspan, delta_t, x0, tolerances, w_k_func = w_func);
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
f = @(t, x, u, p) [x(2); u(1)];
u = @(t, x) -2;
p = 0;

delta_t = 1e-2;

sigma_accel = 0.2; % [m / s2]
G = @(t, x, u, p) [0, 0; 0, sigma_accel] * sqrt(delta_t); % need to double check the sqrt(delta t) part

x0 = [0; 1];

w = @(n) randn([2, n]);

tspan = 0:delta_t:1;

m = 100;

x = zeros([2, numel(tspan), m]);
t = zeros([numel(tspan), m]);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

t_k = tspan;

parfor i = 1:m
    w_k = w(numel(t_k));
    w_func = @(t) w_k(:, floor(t /delta_t) + 1);
    [t(:, i), x(:, :, i)] = sode45(f, G, u, p, w, tspan, delta_t, x0, tolerances, w_k_func = w_func);
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
