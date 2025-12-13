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

delta_t = t_k(2) - t_k(1);

num_traj = size(g_k, 2);

mu = zeros([g_k(1).dim, numel(t_k) - 1, size(twist_k, 3)]);
Y = zeros([g_k(1).dim, g_k(1).dim, numel(t_k) - 1, num_traj]);
for i = 1 : numel(t_k) - 1
    for j = 1 : num_traj
        u_j = u_k_K(:,j) + eps_k(:, i, j);

        mu(:, i, j) = (twist_k(:, i + 1, j) - twist_k(:, i, j)) / delta_t ...
            - f(g_k(i, j), twist_k(:, i, j)) ...
            - B_func(t_k(i), g_k(i, j)) * u_j;

        Y(:, :, i, j) = B_func(t_k(i), g_k(i, j)) * inv(R) * B_func(t_k(i), g_k(i, j))';
    end
end

n = numel(t_k) - 1;
L_k = zeros(n, num_traj);
Ydet_term = zeros(n, num_traj);
for j = 1 : numel(t_k)-1 % loop over time
    for traj_i = 1 : num_traj % loop over trajectories
        u_j = u_k_K(:,j) + eps_k(:, j, traj_i);
        mu_j = mu(:, j, traj_i);
    
        L_k(j, traj_i) = 1/2 * (u_j' * B_func(t_k(j), g_k(j, traj_i))' * inv(Y(:, :, j, traj_i)) * (B_func(t_k(j), g_k(j, traj_i))*u_j + mu_j) * delta_t)...
                        + (1/2 * (mu_j' * inv(Y(:, :, j, traj_i)) * mu_j * delta_t));
        
        Ydet_term(j, traj_i) = det(2 * pi * lambda * Y(:, :, j, traj_i) * delta_t);
    end
end 

L_k = cumsum(L_k, 1, "reverse");
L_k = S_k + L_k + lambda * log(cumprod(sqrt(Ydet_term), 1, "reverse") + 1e-5); % need to add a small constant to avoid log(0)

D = exp(-1 / lambda * (L_k - min(L_k, [], 2)));
D = reshape(D ./ sum(D, 2), 1, n, num_traj);

%D_expec = sum(D .* eps_k / sqrt(delta_t), 3);
D_expec = sum(D .* eps_k, 3);

[~, traj_best] = max(reshape(max(D, [], 2), [], 1));

%u_k_Kp1 = pagemtimes(inv(R) * B' * inv(Y) * B, u_k_K * delta_t + D_expec) / delta_t;
u_k_Kp1 = zeros(size(u_k_K));
for i = 1 : numel(t_k) - 1
    u_k_Kp1(:, i) = inv(R) * B_func(t_k(i), g_k(i, traj_best))' * inv(Y(:, :, i, traj_best)) * B_func(t_k(i), g_k(i, traj_best)) * (u_k_K(:, i) + D_expec(:, i));
end

end