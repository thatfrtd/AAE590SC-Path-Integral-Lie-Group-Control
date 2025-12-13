function [u_new, L_return, L_raw] = iterative_euclidean_update(x_k, v_k, f, B, u_n, delta_u, lambda, R, delta_t, cost_exit, cost_state)

Y = B * inv(R) * B';
Y_uu = R; % B' * inv(Y) * B

L = zeros([size(v_k, 2) - 1, size(v_k, 3)]);
for i = 1 : size(v_k, 3)
    mu_k = diff(v_k(:, :, i), 1, 2) / delta_t - f(:, 2:end, i) - reshape(pagemtimes(B, reshape(u_n + delta_u(:, :, i), 2, 1, [])), 2, []);
    
    % Create Lagrangian
    L(:, i) = cost_exit(i)' ...
        + 1 / 2 * cumsum(cost_state(:, i), 1, "reverse") * delta_t ...
        + 1 / 2 * cumsum(dot(mu_k, inv(Y) * mu_k), 2, "reverse")' * delta_t ...
        + 1 / 2 * cumsum(dot((u_n + delta_u(:, :, i)), (pagemtimes(Y_uu, u_n + delta_u(:, :, i)) + 2 * B' * inv(Y) * mu_k)), 2, "reverse")' * delta_t;
end



% Normalize
w = exp(-(L - min(L, [], 2))/lambda);
Wsum = sum(w, 2);

% Weighted update
dsum = sum(reshape(w, 1, size(w, 1), size(w, 2)) .* delta_u, 3);
u_new = u_n + dsum ./ Wsum';

L_return = w ./ Wsum;
L_raw = L;

end
