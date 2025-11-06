function [t_k, x_k, v_k, u_k, delta_u] = one_step_euler_maruyama_euclidean(f, B, u_k, sigma, t_k, x0, v0)
%ONE_STEP_EULER_MARUYAMA_EUCLIDEAN Summary of this function goes here
%   Detailed explanation goes here
arguments
    f
    B
    u_k
    sigma % pre multiplied by delta t to match SODE45...
    t_k
    x0
    v0
end

w_k = randn(size(u_k));

delta_t = t_k(2) - t_k(1);

% Simulate approximated stochastic differential equation using one step
% Euler Maruyama
x_k = zeros([size(x0, 1), numel(t_k)]);
v_k = zeros([size(v0, 1), numel(t_k)]);
delta_u = zeros([size(u_k, 1), numel(t_k) - 1]);
x_k(:, 1) = x0;
v_k(:, 1) = v0;
for i = 1 : (numel(t_k) - 1)
    x_k(:, i + 1) = x_k(:, i) + v_k(:, i) * delta_t;
    delta_u(:, i) = sigma / sqrt(delta_t) * w_k(:, i);
    v_k(:, i + 1) = v_k(:, i) + f(t_k(i), x_k(:, i)) * delta_t + B(t_k(i), x_k(:, i)) * (u_k(:, i) + delta_u(:, i)) * delta_t;
end
end