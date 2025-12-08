%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 3 December, 2025
% Description: Main script to run the Planar Path Integral Functions on SO3
% Most Recent Change: 4 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = SO3_RotationMatrix();

z = zeros(G.dim);
I = eye(G.dim);

tolerances = odeset(RelTol=1e-6, AbsTol=1e-6);

% Define Noise
noise_delta_t = 1e-2;

w = @(n) randn([G.dim, n]);

sigma_dist = 0.03; % [kg m2 / s2]
sigma = [sigma_dist, 0, 0; 0, sigma_dist, 0; 0, 0, sigma_dist]; % need to double check the sqrt(delta t) part

%% Create Dynamics 
moment_of_inertia = diag([1, 1.11, 1.3]);
J_b = moment_of_inertia; % Generalized inertia

B_map = I;
[f, B] = Euler_Poincare_matrices(SO3_RotationMatrix(), J_b, B_map);

f_with_control = @(t, x, u) [z, I; z, z] * x + [z; B(t, x)] * u;
control_delta_t = 0.01;
t_k = 0:control_delta_t:1.8;
u_k = 0 * ones([G.dim, numel(t_k) - 1]);

R = eye(G.dim);
S_g = 1.6e9 * eye(G.dim);
S_twist = 8e8 * eye(G.dim);

%% Define Initial Condition and Target
g0 = SO3_RotationMatrix(angle2dcm(0, deg2rad(0), deg2rad(0))); % Initial DCM
twist0 = [0; 0; 0]; % Initial angular velocity

gtarg = SO3_RotationMatrix(angle2dcm(0, deg2rad(30), deg2rad(30))'); % Target DCM
twisttarg = zeros([3, 1]);

iterations = 40;
[u_k, cost_t, cost_exit, average_cost] = path_integral_rigidbody_lieSO3(G, z, I, tolerances, noise_delta_t, iterations, sigma, J_b, w, B_map, f, B, f_with_control, control_delta_t, t_k, u_k, R, S_g, S_twist, g0, twist0, gtarg, twisttarg);

% [u_k,  cost_t, cost_exit, average_cost, t_k] = path_integral_spacecraft_lieSO3(G, z, I, tolerances, noise_delta_t, sigma, w, B_map, f, B, f_with_control);

%% Time Histories
sigma_dist = 0.03; % [kg m2 / s2]
sigma = [sigma_dist, 0, 0; 0, sigma_dist, 0; 0, 0, sigma_dist]; % need to double check the sqrt(delta t) part

m = 50;

g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
twist_k = zeros(G.dim, numel(t_k), m); % array of twists (body velocities)
delta_u = zeros([G.dim, numel(t_k) - 1, m]);
w_k = zeros([G.dim, numel(t_k) - 1, m]);

for i = 1 : m
    %[g_k(:, i), twist_k(:, :, i), delta_u(:, :, i), w_k(:,:,i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
    [~, g_k(:, i), twist_k(:, :, i), ~, delta_u(:, :, i)] = sode45_liegroup(u_k, sigma, w, t_k, noise_delta_t, g0, twist0, J_b, tolerances);
end
%[g_k_nom, twist_k_nom, delta_u_nom] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma*0);
[~, g_k_nom, twist_k_nom, ~, delta_u_nom] = sode45_liegroup(u_k, sigma*0, @(n) zeros([G.dim, n]), t_k, noise_delta_t, g0, twist0, J_b, tolerances);
%%
%twist_k_nom = rad2deg(twist_k_nom);
%twist_k = rad2deg(twist_k);

quaternions = SO3_RotationMatrix_to_quaternions(g_k(2:end, :));
quaternions_nom = SO3_RotationMatrix_to_quaternions(g_k_nom);
quaternion_target = SO3_RotationMatrix_to_quaternions(gtarg);

j2 = size(average_cost);

j = j2(2);

figure
tiledlayout(3, 12)

nexttile([1,3])
scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(4, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
scatter(t_k, quaternions_nom(4, :), 10, 'r');
scatter(t_k(end), quaternion_target(4), 10, 'k');
xlabel("Time [s]")
ylabel("q_w")
title("")

nexttile([1,3])
scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(1, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
%         alpha(D_k(:));
scatter(t_k, quaternions_nom(1, :), 10, 'r');
scatter(t_k(end), quaternion_target(1), 10, 'k');
xlabel("Time [s]")
ylabel("q_x")
title("")
nexttile([1,3])
scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(2, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
%         alpha(D_k(:));
scatter(t_k, quaternions_nom(2, :), 10, 'r');
scatter(t_k(end), quaternion_target(2), 10, 'k');
xlabel("Time [s]")
ylabel("q_y")
title("")
nexttile([1,3])
scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(3, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
scatter(t_k, quaternions_nom(3, :), 10, 'r');
scatter(t_k(end), quaternion_target(3), 10, 'k');
xlabel("Time [s]")
ylabel("q_z")
title("")
sgtitle(sprintf("Unit Quaternions: Iteration %d, avg cost %.3f", j, average_cost(j)))

nexttile([1,4])
scatter(repmat(t_k, 1, m), squeeze(twist_k(1, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
%         alpha(D_k(:));
scatter(t_k, twist_k_nom(1, :), 10, 'r');
scatter(t_k(end), twisttarg(1), 10, 'k');
xlabel("Time [s]")
ylabel("w_x")
title("")
nexttile([1,4])
scatter(repmat(t_k, 1, m), squeeze(twist_k(2, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
%         alpha(D_k(:));
scatter(t_k, twist_k_nom(2, :), 10, 'r');
scatter(t_k(end), twisttarg(2), 10, 'k');
xlabel("Time [s]")
ylabel("w_y")
title("")
nexttile([1,4])
scatter(repmat(t_k, 1, m), squeeze(twist_k(3, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
%         alpha(D_k(:));
scatter(t_k, twist_k_nom(3, :), 10, 'r');
scatter(t_k(end), twisttarg(3), 10, 'k');
xlabel("Time [s]")
ylabel("w_x")
title("")

% Control
u_k_m = u_k + delta_u;
nexttile([1,4])
scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(1, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
%         alpha(D_k(:));
scatter(t_k(2:end), u_k(1, :), 10, 'r');
xlabel("Time [s]")
ylabel("u_x")
title("")
nexttile([1,4])
scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(2, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
%         alpha(D_k(:));
scatter(t_k(2:end), u_k(2, :), 10, 'r');
xlabel("Time [s]")
ylabel("u_y")
title("")
nexttile([1,4])
scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(3, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
grid on
%         alpha(D_k(:));
scatter(t_k(2:end), u_k(3, :), 10, 'r');
xlabel("Time [s]")
ylabel("u_z")
title("")


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

%% SO3 Animation
figure
plot_basis(g0.element, "i", ":");
plot_basis(gtarg.element, "f", "--")
clear h
clear g_buffer
clear p_buffer
g_buffer(1) = g0;
g_buffer(2) = g0;
g_buffer(3) = g0;
g_buffer(4) = g0;
g_buffer(5) = g0;
p_buffer(1) = plot_basis(g0.identity, "", "-", scale = 0.1);
h(1) = hgtransform;
set(p_buffer(1).p1,"Parent",h(1))
set(p_buffer(1).p2,"Parent",h(1))
set(p_buffer(1).p3,"Parent",h(1))
p_buffer(2) = plot_basis(g0.identity, "", "-", scale = 0.2);
h(2) = hgtransform;
set(p_buffer(2).p1,"Parent",h(2))
set(p_buffer(2).p2,"Parent",h(2))
set(p_buffer(2).p3,"Parent",h(2))
p_buffer(3) = plot_basis(g0.identity, "", "-", scale = 0.3);
h(3) = hgtransform;
set(p_buffer(3).p1,"Parent",h(3))
set(p_buffer(3).p2,"Parent",h(3))
set(p_buffer(3).p3,"Parent",h(3))
p_buffer(4) = plot_basis(g0.identity, "", "-", scale = 0.4);
h(4) = hgtransform;
set(p_buffer(4).p1,"Parent",h(4))
set(p_buffer(4).p2,"Parent",h(4))
set(p_buffer(4).p3,"Parent",h(4))
p_buffer(5) = plot_basis(g0.identity, "", "-", scale = 1);
h(5) = hgtransform();
set(p_buffer(5).p1,"Parent",h(5))
set(p_buffer(5).p2,"Parent",h(5))
set(p_buffer(5).p3,"Parent",h(5))
for i = 1 : 5 : numel(t_k)
    for j = 1 : 4
        g_buffer(j) = g_buffer(j + 1);
        h(j).Matrix = [g_buffer(j).element, zeros([3, 1]); zeros([1, 3]), 1];
    end
    g_buffer(5) = g_k_nom(i);
    h(5).Matrix = [g_buffer(5).element, zeros([3, 1]); zeros([1, 3]), 1];
    drawnow  % display updates
    pause(control_delta_t / 3)
end
axis equal
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])


%%
figure
plot_basis(g0.element, "i", ":", scale = 1);
for i = 1 : 5 : numel(t_k)
    plot_basis(g_k_nom(i).element, "", "-", scale = 0.7);
end
plot_basis(gtarg.element, "f", "--", scale = 1)
axis equal
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])