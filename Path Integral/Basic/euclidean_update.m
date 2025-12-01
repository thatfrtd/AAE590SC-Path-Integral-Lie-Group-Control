function [u_new, D_k] = euclidean_update(u_n, delta_u, L, lambda)

%normalize
w = exp(-(L - min(L))/lambda);
Wsum = sum(w);

% Weighted update
dsum = sum(reshape(w, 1, 1, []) .* delta_u, 3);
u_new = u_n + dsum / Wsum;

D_k = w ./ Wsum;

end
