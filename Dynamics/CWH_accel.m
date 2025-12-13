function [vdot] = CWH_accel(x, n)
%CWH_ACCEL Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    x
    n
end

r = x(1:3);
v = x(4:6);

vdot = [3 * n ^ 2, 0, 0; 0, 0, 0; 0, 0, -n ^ 2]  * r ...
    +  [0, 2 * n, 0; -2 * n, 0, 0; 0, 0, 0] * v;

end