function [h_func] = chance_constraint_heuristic(c_func, cgrad_func, chess_func, backoff_coef)
%CHANCE_CONSTRAINT_HEURISTIC Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    c_func
    cgrad_func
    chess_func
    backoff_coef
end

h_func = @(x, Ptilde) -c_func(x) - backoff_coef * (sqrt(cgrad_func(x)' * Ptilde * cgrad_func(x)) + sqrt(diag(Ptilde))' * chess_func(x) * sqrt(diag(Ptilde)));

end