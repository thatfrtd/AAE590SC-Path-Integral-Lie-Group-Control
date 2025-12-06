%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 December, 2025
% Description: Main script to run the Planar Path Integral Functions
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

z2 = zeros(2);
I2 = eye(2);

tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

noise_delta_t = 1e-2;

w = @(n) randn([2, n]);

mass = 1;
J_b = mass;

B_map = I2; 
[f, B] = Euler_Poincare_matrices(R2(), J_b, B_map);

f_with_control = @(t, x, u) [z2, I2; z2, z2] * x + [z2; B(t, x)] * u;
control_delta_t = 0.01;
t_k = 0:control_delta_t:1;
u_k = 0 * ones([2, numel(t_k) - 1]);

x0 = [0; 0; 0; 1];
xtarg = [1; 0; 0; 0];
g0 = R2(x0(1:2));
twist0 = x0(3:4);

gtarg = R2(xtarg(1:2));
twisttarg = xtarg(3:4);

u_k = planar_path_integral_linear_lieR2_basic(z2, I2, tolerances, noise_delta_t, w, f, B, f_with_control, control_delta_t, t_k, u_k, x0, xtarg, g0, twist0, gtarg, twisttarg, iterations);

u_k = planar_path_integral_linear_lieR2(z2, I2, tolerances, noise_delta_t, w, f, B, f_with_control, control_delta_t, t_k, u_k, x0, xtarg, g0, twist0, gtarg, twisttarg iterations);