%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 November, 2025
% Description: Test of lie group path integral control on a planar double 
% integrator system represented by the Lie Group R2
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%function u_k = planar_path_integral_linear_lieR2(z2, I2, tolerances, noise_delta_t, w, f, B, f_with_control, control_delta_t, t_k, u_k, x0, xtarg, g0, twist0, gtarg, twisttarg)

G = SE2_RotationMatrix();

z = zeros(G.dim);
I = eye(G.dim);

tolerances = odeset(RelTol=1e-7, AbsTol=1e-7);

%% Define Noise
noise_delta_t = 1e-2;

sigma_torque = 1; % [kg m2 / s2]
sigma_accel = 1; % [kg m / s2]
sigma = diag([sigma_torque*0.001, sigma_accel, sigma_accel]) * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

nu = size(sigma, 1);

w = @(n) randn([nu, n]);

%% Create Dynamics 
mass = 1;
I_b = 1;
J_b = blkdiag(I_b, eye(2) * mass); % Generalized inertia

moment_arm = 1;
B_map = @(g) [1, 0, -moment_arm; zeros([2, 1]), g.R]; 
[f, B] = Euler_Poincare_matrices(G(), J_b, B_map);

f_with_control = @(t, x, u) [z, I; z, z] * x + [z2; B(t, x)] * u;
control_delta_t = 0.01;
t_k = 0:control_delta_t:1;
u_k = 0 * ones([nu, numel(t_k) - 1]);

%% Define Initial Condition and Target
R0 = make_R2(deg2rad(0));
r0 = [0; 0];
w0 = 0;
v0 = [0; 0];
g0 = SE2_RotationMatrix(R0, r0);
twist0 = [w0; v0];

Rtarg = make_R2(deg2rad(180));
rtarg = [1; 0];
wtarg = 0;
vtarg = [0; 0];
gtarg = SE2_RotationMatrix(Rtarg, rtarg);
twisttarg = [wtarg; vtarg];

%% Run Monte Carlo
m = 50;

g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
twist_k = zeros(g0.dim, numel(t_k), m); % array of twists (body velocities)
delta_u = zeros([nu, numel(t_k) - 1, m]);
w_k = zeros([nu, numel(t_k) - 1, m]);

iterations = 150;
R = diag([100, 1, 1])*0.01;
S_g = 2000 * eye(G.dim);
S_twist = 600 * eye(G.dim);
lambda_matrix = sigma * sigma' * R;
lambda = lambda_matrix(1);
average_cost = zeros([1, iterations]);
cost_exit = zeros([m, iterations]);
%%Different Cost_t 
cost_t = zeros([numel(t_k) - 1, m, iterations]);

for j = 1 : iterations
    parfor i = 1:m
        %[t_k2, x_k2, v_k2, u_k2, delta_u2, w_k] = one_step_euler_maruyama_euclidean(@(t,x)0, B, u_k, sigma, t_k, x0(1:2), x0(3:4));
        %[~, x_full, ~, ~] = sode45(f_with_control, u_k, sigma, w, t_k(1:(end - 1)), noise_delta_t, x0, tolerances, w_k = w_k);

        [g_k(:, i), twist_k(:, :, i), delta_u(:, :, i), w_k(:, :, i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);%, w_k = w_k);
    end
    
    % for i = 1:m
    %     [g_k(:, i), twist_k(:, :, i), delta_u(:, :, i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
    % 
    % end

    [L, cost_t(:, :, j), cost_exit(:, j)] = liegroup_cost(g_k, twist_k, u_k + delta_u, control_delta_t, gtarg, twisttarg, R, S_g = S_g, S_twist = S_twist);

%    [u_k, D_k, L_k] = liegroup_update(g_k, u_k, twist_k, pagemtimes(sigma, w_k), f, B, R, L, t_k, lambda);
    [u_k, D_k, L_k] = liegroup_update(g_k, u_k, twist_k, delta_u, f, B, R, L, t_k, lambda);

    average_cost(j) = sum(L_k, "all") / numel(L_k);
    average_cost(j)
    % 
    % if mod(j, 30) == 0 || j == 1
    %     [t_opt, x_opt, ~, delta_u_opt] = sode45(f_with_control, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances);
    % 
    %     figure
    %     x = [g_k(2:end, :).element];
    %     colormap(jet);  
    %     scatter(squeeze(x(1, :)), squeeze(x(2, :)), [], D_k(:), "filled"); hold on
    %     alpha(D_k(:));
    %     plot(x_opt(1, :), x_opt(2, :));
    %     %scatter(xtarg(1), xtarg(2));
    %     xlabel("X [m]")
    %     ylabel("Y [m]")
    %     title("Trajectory")
    %     subtitle(sprintf("Iteration %d, avg cost %.3f", j, average_cost(j)))
    %     grid on
    %     axis equal
    %     colorbar
    % end
end

%% Time Histories
sigma_accel = 0.2; % [m / s2]
sigma = [sigma_accel, 0; 0, sigma_accel] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

m = 50;

g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
twist_k = zeros(G.dim, numel(t_k), m); % array of twists (body velocities)
delta_u = zeros([G.dim, numel(t_k) - 1, m]);
w_k = zeros([G.dim, numel(t_k) - 1, m]);

for i = 1 : m
    [g_k1, twist_k1, delta_u1, w_k] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
    [t, g_k, twist_k, ~, delta_u] = sode45_liegroup(u_k, sigma, w, t_k(1 : (end - 1)), noise_delta_t, g0, twist0, J_b, tolerances, w_k = w_k);
    [t_ck, x_ck, ~, delta_u_ck] = sode45(f_with_control, u_k, sigma, w, t_k(1 : (end - 1)), noise_delta_t, x0, tolerances, w_k = w_k);
end
x = [reshape([g_k.element], G.dim, numel(t_k), m); twist_k];

%[g_k_nom, twist_k_nom, delta_u_nom] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma*0);
[t_nom, g_k_nom, twist_k_nom, u_n, delta_u_nom] = sode45_liegroup(u_k, sigma*0, @(n) zeros([G.dim, n]), t_k, noise_delta_t, g0, twist0, J_b, tolerances);

x_nom = [[g_k_nom.element]; twist_k_nom]; 

% x = zeros([4, numel(t_k), m]);
% t = zeros([numel(t_k), m]);
% delta_u = zeros([2, numel(t_k) - 1, m]);
% 
% for i = 1 : m
%     [t(:, i), x(:, :, i), ~, delta_u(:, :, i)] = sode45(f_with_control, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances);
% end
% [t_nom, x_nom, u_n, ~] = sode45(f_with_control, u_k, sigma*0, w, t_k, noise_delta_t, x0, tolerances);


figure
tiledlayout(1, 3);

nexttile
plot(t, squeeze(x(1, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
plot(t, squeeze(x(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off");
plot(t_nom, squeeze(x_nom(1, :)), Color = "r");
plot(t_nom, squeeze(x_nom(2, :)), Color = "b");
xlabel("Time")
ylabel("Position")
legend("r_x", "r_y")
title("Position vs Time")
grid on

nexttile
plot(t, squeeze(x(3, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
plot(t, squeeze(x(4, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off");
plot(t_nom, squeeze(x_nom(3, :)), Color = "r"); hold on
plot(t_nom, squeeze(x_nom(4, :)), Color = "b");
xlabel("Time")
ylabel("Velocity")
legend("v_x", "v_y")
title("Velocity vs Time")
grid on

nexttile
stairs(t(1 : end - 1, :), squeeze(u_n(1, :) + delta_u(1, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
stairs(t(1 : end - 1, :), squeeze(u_n(2, :) + delta_u(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off");
stairs(t_nom(1 : end - 1, :), squeeze(u_n(1, :)), Color = "r");
stairs(t_nom(1 : end - 1, :), squeeze(u_n(2, :)), Color = "b");
xlabel("Time")
ylabel("Control")
legend("a_x", "a_y")
title("Control Acceleration vs Time")
grid on

%% Trajectory
[g_nom, twist_k_nom] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma * 0);
R_g_nom = createArray([numel(t_k), 1], class(SO2_RotationMatrix));
r_g_nom = zeros([2, numel(t_k)]);
for i = 1 : numel(t_k)
    R_g_nom(i) = SO3_RotationMatrix(g_nom(i).R);
    r_g_nom(:, i) = g_nom(i).r;
end

figure
plot(squeeze(x(1, :, :)), squeeze(x(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility="off"); hold on
plot(squeeze(x_nom(1, :)), squeeze(x_nom(2, :)), Color = "r")
plot(squeeze(x_g_nom(1, :)), squeeze(x_g_nom(2, :)), Color = "g")
xlabel("X [m]")
ylabel("Y [m]")
title("Trajectory")
grid on
axis equal

%% Iterations
figure
plot(average_cost); hold on
scatter(1 : iterations, squeeze(sum(cost_t, 1)), MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
scatter(1 : iterations, cost_exit, MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
plot(mean(squeeze(sum(cost_t, 1)), 1), Color = "b");
plot(mean(cost_exit, 1), Color = "r");
grid on
legend("Total Cost", "Path Cost", "Exit Cost")
xlabel("Iteration")
ylabel("Cost")
title("Cost vs Iteration")
yscale("log")
%end