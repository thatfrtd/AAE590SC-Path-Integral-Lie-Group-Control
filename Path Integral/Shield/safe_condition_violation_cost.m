function [C_safe] = safe_condition_violation_cost(C, h_func, beta, x_k, Ptilde_k)
%SAFE_CONDITION_VIOLATION_COST Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    C
    h_func
    beta
    x_k
    Ptilde_k
end

C_safe = zeros(size(x_k, 2) - 1, size(x_k, 3));
for k = 1 : size(x_k, 2) - 1
    for m = 1 : size(x_k, 3)
        C_safe(k, m) = C * max(-h_func(x_k(:, k, m), Ptilde_k(:, :, k + 1, m)) + (1 - beta) * h_func(x_k(:, k + 1, m), Ptilde_k(:, :, k + 1, m)), 0);
    end
end

end