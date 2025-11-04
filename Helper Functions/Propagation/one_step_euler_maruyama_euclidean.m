function [t, x_k, u_n, delta_u] = one_step_euler_maruyama_euclidean(f, B, u_k, sigma, t_k, x0, tolerances, options)
%ONE_STEP_EULER_MARUYAMA_EUCLIDEAN Summary of this function goes here
%   Detailed explanation goes here
arguments
    f
    B
    u_k
    sigma
    t_k
    x0
    tolerances
    options.w_k_func = []
end

w_k = randn(size(u_k));

delta_t = t_k(2) - t_k(1);

% Simulate approximated stochastic differential equation using one step
% Euler Maruyama
for i = 1 : (numel(t_k) - 1)
    x_k(:, i + 1) = x_k(:, i) + v_k(:, i) * delta_t;
    v_k(:, i + 1) = v_k(:, i) + f(t_k(i), x_k(i)) * delta_t + B(t_k(i), x_k(:, i)) * (u_k(:, i) * delta_t + sigma / sqrt(delta_t) * w_k(:, i));
end
end