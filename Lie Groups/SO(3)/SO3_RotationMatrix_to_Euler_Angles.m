function [euler_angles] = SO3_RotationMatrix_to_Euler_Angles(SO3)
%SO3_ROTATIONMATRIX_TO_EULER_ANGLES Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    SO3
end

N = numel(SO3);
euler_angles = zeros([3, size(SO3)]);

for i = 1 : N
    [euler_angles(1, i), euler_angles(2, i), euler_angles(3, i)] = dcm2angle(SO3(i).element);
end
end