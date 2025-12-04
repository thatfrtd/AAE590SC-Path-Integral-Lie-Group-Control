%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 December, 2025
% Description: Main script to iterating and planar PI eucliden 
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%integration_method = "oneeuler";

z2 = zeros(2);
I2 = eye(2);

A = [z2, I2; z2, z2];
B = [z2; I2];

f = @(t, x, u) pagemtimes(A, x) + pagemtimes(B, u);
control_delta_t = 0.01;
t_k = 0:control_delta_t:1;
u_k = 0 * ones([2, numel(t_k) - 1]);
f_func = @(t, x) z2 * x;
B_func = @(t, x) I2;

x0 = [0; 0];
v0 = [0; 1];
xtarg = [1; 0]; 
vtarg = [0; 0];

noise_delta_t = 1e-2;

w = @(n) randn([2, n]);

tolerances = odeset(RelTol=1e-4, AbsTol=1e-4, InitialStep=0.1, MaxStep=0.1);

u_k = planar_path_integral_linear_iterated_euclidean(integration_method, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, noise_delta_t, w, tolerances);

u_k = planar_path_integral_linear_euclidean(integration_method, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, noise_delta_t, w, tolerances);

u_k = planar_BSS_path_integral_linear_euclidean(integration_method, z2, I2, A, B, f, control_delta_t, t_k, u_k, f_func, B_func, x0, v0, xtarg, vtarg, w, tolerances);