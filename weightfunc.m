function u_new = weightfunc(u_n, delta_u, L, lambda)

%normalize
w = exp(-(L)/lambda);
Wsum = sum(w);

% Weighted update
u_new = u_n;
dsum = sum(w .* delta_u(:, n, :), 2);
u_new(:, n) = u_n(:, n) + (dsum / Wsum).';

end
