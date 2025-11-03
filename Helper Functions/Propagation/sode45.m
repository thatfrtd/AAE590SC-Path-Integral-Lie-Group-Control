function [t, x, u_applied] = sode45(f, G, u, w, tspan, delta_t, x0, tolerances, options)
arguments
    f
    G
    u
    w
    tspan
    delta_t
    x0
    tolerances
    options.w_k_func = []
end
%SODE45 Summary of this function goes here
%   Detailed explanation goes here

if isempty(options.w_k_func)
    w_k = w(numel(tspan(1):delta_t:tspan(end)));
    options.w_k_func = @(t) w_k(:, floor(t /delta_t) + 1);
end

[t, x] = ode45(@(t, x) f(t, x, u(t, x)) + G(t, x, u(t, x)) / sqrt(delta_t) * options.w_k_func(t), tspan, x0, tolerances);

x = x';

u_applied = zeros(size(u(0, x0), 1), numel(t));
for i = 1:numel(t)
    u_applied(:, i) = u(t(i), x(:, i));
end
end

