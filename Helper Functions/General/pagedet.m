function [determinant] = pagedet(A)
%PAGEDET Summary of this function goes here
%   Determinant of a vector of matrices

N = size(A, 3);
determinant = zeros([1, N]);
for i = 1 : N
    determinant(:, :, i) = det(A(:, :, i));
end

end