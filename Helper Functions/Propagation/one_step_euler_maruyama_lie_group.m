function [g_k, twist_k, delta_u, w_k] = one_step_euler_maruyama_lie_group(f, B, g_0, twist_0, u_k, t_k, sigma, options)
%ONE_STEP_EULER_MARUYAMA_LIE_GROUP Summary of this function goes here
%   Detailed explanation goes here
arguments
    f
    B
    g_0
    twist_0
    u_k
    t_k
    sigma % pre multiplied by delta t to match SODE45...
    options.w_k = randn(size(u_k));
end

w_k = options.w_k;

delta_t = t_k(2) - t_k(1);

g_k = createArray(1, numel(t_k), class(g_0)); % array of group elements
g_k(1) = g_0;
twist_k = zeros([g_0.dim, numel(t_k)]); % array of twists (body velocities)
delta_u = zeros([size(u_k, 1), numel(t_k) - 1]);
twist_k(:, 1) = twist_0;
% Simulate approximated stochastic differential equation using one step
% Euler Maruyama on Lie group
for i = 1 : (numel(t_k) - 1)
    g_k(i + 1) = g_k(i).rplus(twist_k(:, i) * delta_t); % reconstruction equation
    delta_u(:, i) = sigma / sqrt(delta_t) * w_k(:, i);
   twist_k(:, i + 1) = twist_k(:, i) + f(g_k(i), twist_k(:, i)) * delta_t + B(t_k(i), g_k(:, i)) * (u_k(:, i) + delta_u(:, i)) * delta_t;
end
end