function [u_k_Kp1] = liegroup_update(g_k, u_k_K, twist_k, eps_k, f, B, R, S_k, t_k, lambda)
%LIEGROUP_UPDATE Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    g_k
    u_k_K
    twist_k
    eps_k
    f
    B
    R
    S_k
    t_k
    lambda
end

delta_t = t_k(2) - t_k(1);

i = 1 : (numel(t_k) - 1); % control indices

Y = B * R \ B';
mu = diff(twist_k, 1, 2) / delta_t ...
    - f(t_k(i), [g_k(i).element]) ...
    - B * u_k_K;

L_k = S_k + 1 / 2 * sum(u_k_K(:, i)' * B' * pagemldivide(Y, (B * u_k_K(:, i) + mu) * delta_t), 3) ...
    + 1 / 2 * sum(pagemtimes(pagetranspose(mu), pagemldivide(Y, mu * delta_t)), 3) ...
    + lambda * log(prod(pagedet(2 * pi * lambda * Y * delta_t))); % Isn't this term constant

D = exp(-1 / lamda * L_k);
D = D / sum(D);

D_expec = sum(D .* eps_k * sqrt(delta_t));

u_k_Kp1 = R \ B' * Y \ B * (u_k_K * delta_t + D_expec);
end