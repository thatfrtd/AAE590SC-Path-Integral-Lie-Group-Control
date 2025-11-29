function u_new = euclidean_update(u_n, delta_u, L, lambda)

%normalize
w = exp(-(L - min(L))/lambda);
Wsum = sum(w);

% Weighted update
dsum = sum(reshape(w, 1, 1, []) .* delta_u, 3);
u_new = u_n + dsum / Wsum;

end
