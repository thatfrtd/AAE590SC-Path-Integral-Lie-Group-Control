function [u_k_Kp1, D, L_k] = liegroup_update(g_k, u_k_K, twist_k, eps_k, f, B_func, R, S_k, t_k, lambda)
%LIEGROUP_UPDATE Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    g_k
    u_k_K
    twist_k
    eps_k
    f
    B_func
    R
    S_k
    t_k
    lambda
end

B = B_func([], []); % ASSUMING B IS CONSTANT

delta_t = t_k(2) - t_k(1);

Y = B * R \ B';
mu = zeros([g_k(1).dim, numel(t_k) - 1, size(twist_k, 3)]);

for i = 1 : numel(t_k) - 1
    mu(:, i, :) = (twist_k(:, i + 1, :) - twist_k(:, i, :)) / delta_t ...
        - f(g_k(i, :), twist_k(:, i, :)) ...
        - B * u_k_K(:, i);
end

% L_k = S_k + 1 / 2 * sum(u_k_K(:, i)' * B' * pagemldivide(Y, (B * u_k_K(:, i) + mu) * delta_t), 3) ...
%     + 1 / 2 * sum(pagemtimes(pagetranspose(mu), pagemldivide(Y, mu * delta_t)), 3) ...
%     + lambda * log(prod(sqrt(pagedet(2 * pi * lambda * Y * delta_t)))); % Isn't this term constant
n = numel(t_k) - 1;
num_traj = size(g_k, 2);
L_k = zeros(n, num_traj);
for j = 1 : numel(t_k)-1 % loop over time
    for traj_i = 1 : num_traj % loop over trajectories
        u_j = u_k_K(:,j);
        mu_j = mu(:, j, traj_i);
    
        L_k(j, traj_i) = 1/2 * (u_j' * B' * inv(Y) * (B*u_j + mu_j) * delta_t)...
                        + (1/2 * (mu_j' * inv(Y) * mu_j * delta_t));
    end
end 

L_k = cumsum(L_k, 1, "reverse");
L_k = S_k + L_k*0 + 0*lambda * log(cumprod(sqrt(det(2 * pi * lambda * Y * delta_t))));

D = exp(-1 / lambda * (L_k - min(L_k, [], 2)));
%D = sum(reshape(D ./ sum(D, 2), 1, n, num_traj), 2) / size(D, 1);
D = reshape(D ./ sum(D, 2), 1, n, num_traj);% / size(D, 1);

D_expec = sum(D .* eps_k * sqrt(delta_t), 3);

u_k_Kp1 = pagemtimes(inv(R) * B' * inv(Y) * B, u_k_K * delta_t + D_expec) / delta_t;

end