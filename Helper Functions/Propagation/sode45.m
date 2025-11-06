function [t, x, u_n, delta_u] = sode45(f, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances, options)
arguments
    f
    u_k
    sigma
    w
    t_k
    noise_delta_t
    x0
    tolerances
    options.w_k_func = []
end
%SODE45 Summary of this function goes here
%   Detailed explanation goes here

% Create continuous noise function
if isempty(options.w_k_func)
    w_k = w(numel(t_k(1):noise_delta_t:t_k(end)));
    options.w_k_func = @(t) w_k(:, floor(t /noise_delta_t) + 1);
end

% Create continuous control function (assume zero order hold)
control_delta_t = t_k(2) - t_k(1);
u_k_func = @(t) u_k(:, min(floor(t / control_delta_t) + 1, size(u_k, 2)));

% Simulate approximated stochastic differential equation
[t, x] = ode45(@(t, x) f(t, x, u_k_func(t) + sigma / sqrt(noise_delta_t) * options.w_k_func(t)), t_k, x0, tolerances);

% I like state as first dimension and timestep as second
x = x';

% Get continuous control array
u_n = zeros(size(u_k, 1), numel(t) - 1);
delta_u = zeros(size(u_k, 1), numel(t) - 1);
for i = 1:(numel(t) - 1)
    u_n(:, i) = u_k_func(t(i));
    delta_u(:, i) = sigma / sqrt(noise_delta_t) * options.w_k_func(t(i));
end
end

