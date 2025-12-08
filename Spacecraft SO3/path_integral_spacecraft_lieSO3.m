%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 22 November, 2025
% Description: Lie group path integral control of a spacecraft's
% 3D orientation which is in a circular orbit with a gravity gradient
% moment
% Most Recent Change: 22 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [u_k, cost_t, cost_exit, average_cost, t_k] = path_integral_spacecraft_lieSO3(G, z, I, tolerances, noise_delta_t,sigma, w, B_map, f, B, f_with_control)

%G = SO3_RotationMatrix();

% z = zeros(G.dim);
% I = eye(G.dim);
% 
% tolerances = odeset(RelTol=1e-12, AbsTol=1e-12);

%% Define Noise
%noise_delta_t = 2e-1;

%w = @(n) randn([G.dim, n]);

% sigma_dist = 15000; % [kg m2 / s2]
% sigma = [sigma_dist/10, 0, 0; 0, sigma_dist, 0; 0, 0, sigma_dist] * sqrt(noise_delta_t); % need to double check the sqrt(delta t) part

%% Create Dynamics 
% moment_of_inertia = [1.29 0 0; 0 9.68 0; 0 0 10.1] * 1e6; % [kg m2]
% J_b = moment_of_inertia; % Generalized inertia

% Create orbit dynamics
R_E = 6378.1363; % [km] 
h = 300; % [km]
mu_E = 398600.4415; % [km3 / s2]
mean_motion = sqrt(mu_E / (R_E + h) ^ 3); % [rad / s]
M_grav_gradient = @(x) gravity_gradient(x, mean_motion, moment_of_inertia);

% Create Euler Poincare dynamics
B_map = I*1e-0;
[f, B] = Euler_Poincare_matrices(SO3_RotationMatrix(), J_b, B_map, M_grav_gradient);

f_with_control = @(t, x, u) [z, I; z, z] * x + [z; B(t, x)] * u;
control_delta_t = 0.2;
t_k = 0:control_delta_t:50;
u_k = 0 * ones([G.dim, numel(t_k) - 1]);

R = 1e-3*eye(G.dim);
S_g = 1.7e6 * eye(G.dim);
S_twist = 8e5 * eye(G.dim);


%% Define Initial Condition and Ta rget
g0 = SO3_RotationMatrix(angle2dcm(0.1, 1, 0.4)); % Initial DCM
twist0 = [0; 0.1; 0]; % Initial angular velocity
xtarg = G.identity;

gtarg = SO3_RotationMatrix(eye(3)); % Identity element
twisttarg = zeros([3, 1]);

%% Run Monte Carlo
m = 50;

g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
twist_k = zeros(G.dim, numel(t_k), m); % array of twists (body velocities)
delta_u = zeros([G.dim, numel(t_k) - 1, m]);
w_k = zeros([G.dim, numel(t_k) - 1, m]);


iterations = 100;
lambda_matrix = sigma * sigma' * R;
lambda = lambda_matrix(1);
average_cost = zeros([1, iterations]);
cost_t = zeros([numel(t_k) - 1, m, iterations]);
cost_exit = zeros([m, iterations]);
%%
for j = 1 : iterations
    parfor i = 1 : m
        [g_k(:, i), twist_k(:, :, i), delta_u(:, :, i), w_k(:,:,i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
    end
    
    [L, cost_t(:, :, j), cost_exit(:, j)] = liegroup_cost(g_k, twist_k, u_k + delta_u, control_delta_t, gtarg, twisttarg, R*0, S_g = S_g, S_twist = S_twist);
    [u_k, D_k] = liegroup_update(g_k, u_k, twist_k, delta_u, f, B, R, L, t_k, lambda);

    average_cost(j) = sum(L .* control_delta_t, "all") / size(L, 2);
    average_cost(j)

    if mod(j, 20) == 0 || j == 1
        %[t_opt, x_opt, ~, delta_u_opt] = sode45(f_with_control, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances);

        figure
        tiledlayout(3, 12)
        
        quaternions = SO3_RotationMatrix_to_quaternions(g_k(2:end, :));
        quaternion_target = SO3_RotationMatrix_to_quaternions(gtarg);
        colormap(jet);  
        
        nexttile([1,3])
        scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(4, :)), [], L(:), "filled"); hold on
        grid on
         alpha(D_k(:));
         scatter(t_k(end), quaternion_target(4), 10, 'r');
        colorbar
        xlabel("Time [s]")
        ylabel("q_w")
        title("")
        nexttile([1,3])
        scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(1, :)), [], L(:), "filled"); hold on
        grid on
         alpha(D_k(:));
        scatter(t_k(end), quaternion_target(1), 10, 'r');
        colorbar
        xlabel("Time [s]")
        ylabel("q_x")
        title("")
        nexttile([1,3])
        scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(2, :)), [], L(:), "filled"); hold on
        grid on
         alpha(D_k(:));
        scatter(t_k(end), quaternion_target(2), 10, 'r');
        colorbar
        xlabel("Time [s]")
        ylabel("q_y")
        title("")
        nexttile([1,3])
        scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(3, :)), [], L(:), "filled"); hold on
        grid on
         alpha(D_k(:));
        scatter(t_k(end), quaternion_target(3), 10, 'r');
        colorbar
        xlabel("Time [s]")
        ylabel("q_z")
        title("")
        sgtitle(sprintf("Unit Quaternions: Iteration %d, avg cost %.3f", j, average_cost(j)))

        twist_k_Nu = twist_k(:, 2:end, :);

        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(twist_k_Nu(1, :)), [], L(:), "filled"); hold on
%         alpha(D_k(:));
        grid on
        colorbar
        xlabel("Time [s]")
        ylabel("w_x")
        title("")
        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(twist_k_Nu(2, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
        colorbar
        xlabel("Time [s]")
        ylabel("w_y")
        title("")
        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(twist_k_Nu(3, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
        colorbar
        xlabel("Time [s]")
        ylabel("w_z")
        title("")

        % Control
        u_k_m = u_k + delta_u;

        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(1, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
        colorbar
        xlabel("Time [s]")
        ylabel("u_x")
        title("")
        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(2, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
        colorbar
        xlabel("Time [s]")
        ylabel("u_y")
        title("")
        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(3, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
        colorbar
        xlabel("Time [s]")
        ylabel("u_z")
        title("")
    end
end


% %% Time Histories
% sigma_dist = 0.001; % [kg m2 / s2]
% sigma = [sigma_dist, 0, 0; 0, sigma_dist, 0; 0, 0, sigma_dist]; % need to double check the sqrt(delta t) part
% 
% m = 50;
% 
% g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
% twist_k = zeros(G.dim, numel(t_k), m); % array of twists (body velocities)
% delta_u = zeros([G.dim, numel(t_k) - 1, m]);
% w_k = zeros([G.dim, numel(t_k) - 1, m]);
% 
% for i = 1 : m
%     [g_k(:, i), twist_k(:, :, i), delta_u(:, :, i), w_k(:,:,i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
% end
% [g_k_nom, twist_k_nom, delta_u_nom, w_k_nom] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma * 0);
% twist_k_nom = rad2deg(twist_k_nom);
% twist_k = rad2deg(twist_k);
% 
% quaternions = SO3_RotationMatrix_to_quaternions(g_k(2:end, :));
% quaternions_nom = SO3_RotationMatrix_to_quaternions(g_k_nom);
% quaternion_target = SO3_RotationMatrix_to_quaternions(gtarg);
% 
% 
% figure
% tiledlayout(3, 12)
% 
% nexttile([1,3])
% scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(4, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% scatter(t_k, quaternions_nom(4, :), 10, 'r');
% scatter(t_k(end), quaternion_target(4), 10, 'k');
% xlabel("Time [s]")
% ylabel("q_w")
% title("")
% 
% nexttile([1,3])
% scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(1, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% %         alpha(D_k(:));
% scatter(t_k, quaternions_nom(1, :), 10, 'r');
% scatter(t_k(end), quaternion_target(1), 10, 'k');
% xlabel("Time [s]")
% ylabel("q_x")
% title("")
% nexttile([1,3])
% scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(2, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% %         alpha(D_k(:));
% scatter(t_k, quaternions_nom(2, :), 10, 'r');
% scatter(t_k(end), quaternion_target(2), 10, 'k');
% xlabel("Time [s]")
% ylabel("q_y")
% title("")
% nexttile([1,3])
% scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(3, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% scatter(t_k, quaternions_nom(3, :), 10, 'r');
% scatter(t_k(end), quaternion_target(3), 10, 'k');
% xlabel("Time [s]")
% ylabel("q_z")
% title("")
% sgtitle(sprintf("Unit Quaternions: Iteration %d, avg cost %.3f", j, average_cost(j)))
% 
% nexttile([1,4])
% scatter(repmat(t_k, 1, m), squeeze(twist_k(1, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% %         alpha(D_k(:));
% scatter(t_k, twist_k_nom(1, :), 10, 'r');
% scatter(t_k(end), twisttarg(1), 10, 'k');
% xlabel("Time [s]")
% ylabel("w_x")
% title("")
% nexttile([1,4])
% scatter(repmat(t_k, 1, m), squeeze(twist_k(2, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% %         alpha(D_k(:));
% scatter(t_k, twist_k_nom(2, :), 10, 'r');
% scatter(t_k(end), twisttarg(2), 10, 'k');
% xlabel("Time [s]")
% ylabel("w_y")
% title("")
% nexttile([1,4])
% scatter(repmat(t_k, 1, m), squeeze(twist_k(3, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% %         alpha(D_k(:));
% scatter(t_k, twist_k_nom(3, :), 10, 'r');
% scatter(t_k(end), twisttarg(3), 10, 'k');
% xlabel("Time [s]")
% ylabel("w_x")
% title("")
% 
% % Control
% u_k_m = u_k + delta_u;
% nexttile([1,4])
% scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(1, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% %         alpha(D_k(:));
% scatter(t_k(2:end), u_k(1, :), 10, 'r');
% xlabel("Time [s]")
% ylabel("u_x")
% title("")
% nexttile([1,4])
% scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(2, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% %         alpha(D_k(:));
% scatter(t_k(2:end), u_k(2, :), 10, 'r');
% xlabel("Time [s]")
% ylabel("u_y")
% title("")
% nexttile([1,4])
% scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(3, :)), [], [192, 192, 192] / 256, "filled", HandleVisibility="off"); hold on
% grid on
% %         alpha(D_k(:));
% scatter(t_k(2:end), u_k(3, :), 10, 'r');
% xlabel("Time [s]")
% ylabel("u_z")
% title("")
% 
% 
% 
% %% Iterations
% figure
% plot(average_cost); hold on
% scatter(1 : iterations, squeeze(sum(cost_t, 1)), MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
% scatter(1 : iterations, cost_exit, MarkerFaceColor = [192, 192, 192] / 256, MarkerEdgeColor="none", HandleVisibility="off")
% plot(mean(squeeze(sum(cost_t, 1)), 1), Color = "b");
% plot(mean(cost_exit, 1), Color = "r");
% grid on
% legend("Total Cost", "Path Cost", "Exit Cost")
% xlabel("Iteration")
% ylabel("Cost")
% title("Cost vs Iteration")
% yscale("log")
% 
% %% SO3 Animation
% figure
% plot_basis(g0.element, "i", ":");
% plot_basis(gtarg.element, "f", "--")
% clear h
% clear g_buffer
% clear p_buffer
% g_buffer(1) = g0;
% g_buffer(2) = g0;
% g_buffer(3) = g0;
% g_buffer(4) = g0;
% g_buffer(5) = g0;
% p_buffer(1) = plot_basis(g0.identity, "", "-", scale = 0.1);
% h(1) = hgtransform;
% set(p_buffer(1).p1,"Parent",h(1))
% set(p_buffer(1).p2,"Parent",h(1))
% set(p_buffer(1).p3,"Parent",h(1))
% p_buffer(2) = plot_basis(g0.identity, "", "-", scale = 0.2);
% h(2) = hgtransform;
% set(p_buffer(2).p1,"Parent",h(2))
% set(p_buffer(2).p2,"Parent",h(2))
% set(p_buffer(2).p3,"Parent",h(2))
% p_buffer(3) = plot_basis(g0.identity, "", "-", scale = 0.3);
% h(3) = hgtransform;
% set(p_buffer(3).p1,"Parent",h(3))
% set(p_buffer(3).p2,"Parent",h(3))
% set(p_buffer(3).p3,"Parent",h(3))
% p_buffer(4) = plot_basis(g0.identity, "", "-", scale = 0.4);
% h(4) = hgtransform;
% set(p_buffer(4).p1,"Parent",h(4))
% set(p_buffer(4).p2,"Parent",h(4))
% set(p_buffer(4).p3,"Parent",h(4))
% p_buffer(5) = plot_basis(g0.identity, "", "-", scale = 1);
% h(5) = hgtransform();
% set(p_buffer(5).p1,"Parent",h(5))
% set(p_buffer(5).p2,"Parent",h(5))
% set(p_buffer(5).p3,"Parent",h(5))
% for i = 1 : 10 : numel(t_k)
%     for j = 1 : 4
%         g_buffer(j) = g_buffer(j + 1);
%         h(j).Matrix = [g_buffer(j).element, zeros([3, 1]); zeros([1, 3]), 1];
%     end
%     g_buffer(5) = g_k_nom(i);
%     h(5).Matrix = [g_buffer(5).element, zeros([3, 1]); zeros([1, 3]), 1];
%     drawnow  % display updates
%     pause(control_delta_t / 3)
% end
% axis equal
% xlim([-1, 1])
% ylim([-1, 1])
% zlim([-1, 1])