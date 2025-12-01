function [quaternions] = SO3_RotationMatrix_to_quaternions(SO3)
%SO3_ROTATIONMATRIX_TO_EULER_ANGLES Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    SO3
end

N = numel(SO3);
quaternions = zeros([4, size(SO3)]);

for i = 1 : N
    quaternions(:, i) = SO3_quaternion.Exp(SO3(i).Log()).element;
end
end