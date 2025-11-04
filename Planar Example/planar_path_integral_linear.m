%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Travis Hastreiter 
% Created On: 3 November, 2025
% Description: Test of stochastic propagation without optimized control
% Most Recent Change: 3 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z2 = zeros(2);
I2 = eye(2);

f = @(t, x, u) A * x + B * (u);
t_k = 0:0.01:1;
u_k = -2 * ones([2, numel(t_k)]);

noise_delta_t = 1e-2;

w = @(n) randn([2, n]);

sigma_accel = 0.2; % [m / s2]
G = @(t, x, u) [0, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

x0 = [0; 1];

m = 100;

x = zeros([2, numel(tspan), m]);
t = zeros([numel(tspan), m]);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

parfor i = 1:m
    [t(:, i), x(:, :, i)] = sode45(f, G, u_k, w, t_k, noise_delta_t, x0, tolerances);
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
