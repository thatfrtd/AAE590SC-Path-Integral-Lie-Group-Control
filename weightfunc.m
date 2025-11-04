function u_new = weightfunc(u_n, deltaU, L, lambda)

%normalize
w = exp(-(L)/lambda);
Wsum = sum(w);

% Weighted update
u_new = u_n;
for n = 1:N
    dsum = zeros(1, Nu);
    for k = 1:K
        dsum = dsum + w(k) * (deltaU(k,n,:)).';
    end
    u_new(:, n) = u_n(:, n) + (dsum / Wsum).';
end

end
