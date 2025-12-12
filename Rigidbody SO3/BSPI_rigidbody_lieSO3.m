%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% AAE 590SC
% Path Integral Control on Lie Groups
% Author: Nyssa Guha, Travis Hastreiter
% Created On: 21 November, 2025
% Description: Lie group path integral control of a rigid body's
% 3D orientation
% Most Recent Change: 21 November, 2025
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

G = SO3_RotationMatrix();

z = zeros(G.dim);
I = eye(G.dim);

tolerances = odeset(RelTol=1e-6, AbsTol=1e-6);

%% Define Noise
noise_delta_t = 1e-2;

w = @(n) randn([G.dim, n]);

sigma_dist = 0.03; % [kg m2 / s2]
sigma = [sigma_dist, 0, 0; 0, sigma_dist, 0; 0, 0, sigma_dist]; % need to double check the sqrt(delta t) part

%% Create Dynamics 
moment_of_inertia = diag([1, 1.11, 1.3]);
J_b = moment_of_inertia; % Generalized inertia

B_map = @(x) I;
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

gtarg = SO3_RotationMatrix(angle2dcm(50, deg2rad(120), deg2rad(260))'); % Target DCM
twisttarg = zeros([3, 1]);

%%
% Define initial state covariance
sigma_xhat0 = [5e-2; ... % _x
               5e-2; ... % _y
               5e-2]; ... % _z
P0 = diag(sigma_xhat0 .^ 2)*1e-5;

P_k = zeros([3, 3, numel(t_k)]);
P_k(:, :, 1) = P0;
X0 = chol(P0);

filter.A_k = G.Ad(G.Exp(twist0 * control_delta_t));
filter.G_k = B(0, g0) * sigma / sqrt(noise_delta_t) * control_delta_t;
for k = 1 : numel(t_k) - 1
    P_k(:, :, k + 1) = filter.A_k * P_k(:, :, k) * filter.A_k' ...
                          + filter.G_k * filter.G_k';
end

P_fail = 0.01;
backoff_coef = sqrt(chi2inv(1 - P_fail, G.dim));

% %%
% g0.sample_left_gaussian()
% 
% figure
% plot_basis(g0.element, "i", ":", scale = 1);
% for i = 1 : 5 : numel(t_k)
%     plot_basis(g_k_nom(i).element, "", "-", scale = 0.7);
% end
% plot_basis(gtarg.element, "f", "--", scale = 1)
% axis equal
% xlim([-1, 1])
% ylim([-1, 1])
% zlim([-1, 1])

%% Set up constraint
logistic = @(x) 1 ./ (1 + exp(-8000*x)); % Smooth out penalty

objective_penalty = 1e6;
cone_angle = 30; % [deg]
sun_dcm = angle2dcm(-1.4, 0.8, 1.2);
sun_direction = sun_dcm * [0; 0; 1];
%sensor_direction_body = angle2dcm(-1, 0.2, -2) * [0; 0; 1];
sensor_direction_body = angle2dcm(0, 0, 0) * [0; 1; 0];
sensor_direction_inertial = gtarg.element * sensor_direction_body;

alpha = @(g, X) X' * g.element' * [-(sun_direction(2) * sensor_direction_inertial(3) - sun_direction(3) * sensor_direction_inertial(2)); ...
                                     sun_direction(1) * sensor_direction_inertial(3) - sun_direction(3) * sensor_direction_inertial(1); ...
                                   -(sun_direction(1) * sensor_direction_inertial(2) - sun_direction(2) * sensor_direction_inertial(1))];
alpha_hat = @(g, X) alpha(g, X) / norm(alpha(g, X));
sun_chance_constraint = @(g, X) -cosd(cone_angle) + sun_direction' * g.rplus(backoff_coef * X * alpha_hat(g, X) * 0).element * sensor_direction_body;
chance_constrained_state_running_cost = @(g, twist, X) objective_penalty * logistic(sun_chance_constraint(g, X));

sun_constraint = @(g, twist) objective_penalty * logistic(-cosd(cone_angle) + sun_direction' * g.element * sensor_direction_body);


%% Run Monte Carlo
m = 40;

g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
twist_k = zeros(G.dim, numel(t_k), m); % array of twists (body velocities)
P_k = zeros(G.dim, G.dim, numel(t_k), m);
P_k(:, :, 1, :) = repmat(P0, 1, 1, 1, m);
X_k = zeros(G.dim, G.dim, numel(t_k), m);
X_k(:, :, 1, :) = repmat(chol(P0), 1, 1, 1, m);
delta_u = zeros([G.dim, numel(t_k) - 1, m]);
w_k = zeros([G.dim, numel(t_k) - 1, m]);


iterations = 120;
lambda_matrix = sigma * sigma' * R;
lambda = lambda_matrix(1) * 100;
average_cost = zeros([1, iterations]);
cost_t = zeros([numel(t_k) - 1, m, iterations]);
cost_exit = zeros([m, iterations]);
for j = 1 : iterations
    %[g_k(:, 1), twist_k(:, :, 1), delta_u(:, :, 1), w_k(:,:,1)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma*0);
    %[~, g_k(:, 1), twist_k(:, :, 1), ~, delta_u(:, :, 1)] = sode45_liegroup(u_k, sigma*0, w, t_k, noise_delta_t, g0, twist0, J_b, tolerances);
    for i = 1 : m
        P_k_i = zeros([3, 3, numel(t_k)]);
        P_k_i(:, :, 1) = P0;
        X_k_i = zeros([3, 3, numel(t_k)]);
        X_k_i(:, :, 1) = X0;
        %[~, g_k(:, i), twist_k(:, :, i), ~, delta_u(:, :, i)] = sode45_liegroup(u_k, sigma, w, t_k, noise_delta_t, g0, twist0, J_b, tolerances);
        [g_k_i, twist_k_i, delta_u(:, :, i), w_k(:,:,i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
        for k = 1 : numel(t_k) - 1
            A_k = G.Ad(G.Exp(twist_k_i(:, k) * control_delta_t));
            G_k = B(t_k(k), g_k_i(k)) * sigma / sqrt(noise_delta_t) * control_delta_t;
            
            P_k_i(:, :, k + 1) = A_k * P_k_i(:, :, k) * A_k' + G_k * G_k';        
            X_k_i(:, :, k + 1) = chol(P_k_i(:, :, k + 1));
        end

        g_k(:, i) = g_k_i;
        twist_k(:, :, i) = twist_k_i;
        P_k(:, :, :, i) = P_k_i;
        X_k(:, :, :, i) = X_k_i;
    end
    
    [L, cost_t(:, :, j), cost_exit(:, j)] = liegroup_cost(g_k, twist_k, u_k + delta_u, control_delta_t, gtarg, twisttarg, R, S_g = S_g, S_twist = S_twist);
    %[L, cost_t(:, :, j), cost_exit(:, j)] = liegroup_cost(g_k, twist_k, u_k + delta_u, control_delta_t, gtarg, twisttarg, R, S_g = S_g, S_twist = S_twist, state_running_cost = sun_constraint);

%    [u_k, D_k] = liegroup_update(g_k, u_k, twist_k, pagemtimes(sigma, w_k), f, B, R, L, t_k, lambda);
    u_km1 = u_k;
    S_k = cumsum(L, 1, "reverse");
    [u_k, D_k] = liegroup_update(g_k, u_k, twist_k, delta_u, f, B, R, S_k, t_k, lambda);

    average_cost(j) = min(sum(L .* control_delta_t, 1), [], "all");
    average_cost(j)

    if mod(j, 40) == 0 || j == 1
        %[t_opt, x_opt, ~, delta_u_opt] = sode45(f_with_control, u_k, sigma, w, t_k, noise_delta_t, x0, tolerances);

        figure
        tiledlayout(3, 12)
        
        quaternions = SO3_RotationMatrix_to_quaternions(g_k(2:end, :));
        quaternion_target = SO3_RotationMatrix_to_quaternions(gtarg);
        colormap(jet);  
        
        nexttile([1,3])
        scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(4, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
         scatter(t_k(end), quaternion_target(4), 10, 'r');
        colorbar
        xlabel("Time [s]")
        ylabel("q_w")
        title("")
        nexttile([1,3])
        scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(1, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
        scatter(t_k(end), quaternion_target(1), 10, 'r');
        colorbar
        xlabel("Time [s]")
        ylabel("q_x")
        title("")
        nexttile([1,3])
        scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(2, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
        scatter(t_k(end), quaternion_target(2), 10, 'r');
        colorbar
        xlabel("Time [s]")
        ylabel("q_y")
        title("")
        nexttile([1,3])
        scatter(repmat(t_k(2:end), 1, m), squeeze(quaternions(3, :)), [], L(:), "filled"); hold on
        grid on
%         alpha(D_k(:));
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
        u_k_m = u_km1 + delta_u;

        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(1, :)), [], L(:), "filled"); hold on
        scatter(t_k(2:end), squeeze(u_k(1, :)), 50, "black", "square"); hold on
        grid on
%         alpha(D_k(:));
        colorbar
        xlabel("Time [s]")
        ylabel("u_x")
        title("")
        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(2, :)), [], L(:), "filled"); hold on
        scatter(t_k(2:end), squeeze(u_k(2, :)), 50, "black", "square"); hold on
        grid on
%         alpha(D_k(:));
        colorbar
        xlabel("Time [s]")
        ylabel("u_y")
        title("")
        nexttile([1,4])
        scatter(repmat(t_k(2:end), 1, m), squeeze(u_k_m(3, :)), [], L(:), "filled"); hold on
        scatter(t_k(2:end), squeeze(u_k(3, :)), 50, "black", "square"); hold on
        grid on
%         alpha(D_k(:));
        colorbar
        xlabel("Time [s]")
        ylabel("u_z")
        title("")
    end
end


%% Time Histories
sigma_dist = 0.03; % [kg m2 / s2]
sigma = [sigma_dist, 0, 0; 0, sigma_dist, 0; 0, 0, sigma_dist]; % need to double check the sqrt(delta t) part

m = 50;

g_k = createArray(numel(t_k), m, class(g0)); % array of group elements
twist_k = zeros(G.dim, numel(t_k), m); % array of twists (body velocities)
delta_u = zeros([G.dim, numel(t_k) - 1, m]);
w_k = zeros([G.dim, numel(t_k) - 1, m]);

parfor i = 1 : m
    %[g_k(:, i), twist_k(:, :, i), delta_u(:, :, i), w_k(:,:,i)] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma);
    [~, g_k(:, i), twist_k(:, :, i), ~, delta_u(:, :, i)] = sode45_liegroup(u_k, sigma, w, t_k, noise_delta_t, g0, twist0, J_b, tolerances);
end
%[g_k_nom, twist_k_nom, delta_u_nom] = one_step_euler_maruyama_lie_group(f, B, g0, twist0, u_k, t_k, sigma*0);
[~, g_k_nom, twist_k_nom, ~, delta_u_nom] = sode45_liegroup(u_k, sigma*0, @(n) zeros([G.dim, n]), t_k, noise_delta_t, g0, twist0, J_b, tolerances);

%%
lerror = zeros([3, numel(t_k), m]);
for i = 1 : m
    for k = 1 : numel(t_k)
        lerror(:, k, i) = g_k_nom(k).left_invariant_error(g_k(k, i)); 
    end
end

for k = 1 : numel(t_k) - 1
    filter.A_k = G.Ad(G.Exp(twist_k_nom(:, k) * control_delta_t));
    filter.G_k = B(t_k(k), g_k_nom(k)) * sigma / sqrt(noise_delta_t) * control_delta_t;
    
    P_k(:, :, k + 1) = filter.A_k * P_k(:, :, k) * filter.A_k' ...
                          + filter.G_k * filter.G_k';
end

error_std = zeros([3, numel(t_k)]);
for k = 1 : numel(t_k)
    error_std(:, k) = rad2deg(diag(chol(P_k(:, :, k)))) * backoff_coef;
end
%%

tiledlayout(1, 3)
nexttile
plot(t_k, rad2deg((squeeze(lerror(1, :, :)))), Color = [192, 192, 192] / 256); hold on 
plot(t_k, error_std(1, :), Color = "m")
plot(t_k, -error_std(1, :), Color = "m")
grid on
nexttile
plot(t_k, rad2deg((squeeze(lerror(2, :, :)))), Color = [192, 192, 192] / 256); hold on
plot(t_k, error_std(2, :), Color = "m")
plot(t_k, -error_std(2, :), Color = "m")
grid on
nexttile
plot(t_k, rad2deg((squeeze(lerror(3, :, :)))), Color = [192, 192, 192] / 256); hold on
plot(t_k, error_std(3, :), Color = "m")
plot(t_k, -error_std(3, :), Color = "m")
grid on

%%
%twist_k_nom = rad2deg(twist_k_nom);
%twist_k = rad2deg(twist_k);

quaternions = SO3_RotationMatrix_to_quaternions(g_k(2:end, :));
quaternions_nom = SO3_RotationMatrix_to_quaternions(g_k_nom);
quaternion_target = SO3_RotationMatrix_to_quaternions(gtarg);


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
% for i = 1 : 5 : numel(t_k)
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
% 

%%
figure
plot_basis(g0.element, "i", ":", scale = 1);
for i = 1 : 5 : numel(t_k)
    plot_basis(g_k_nom(i).element, "", "-", scale = 0.7);
    sensor_direction_inertial = g_k_nom(i).element * sensor_direction_body;
    quiver3(0, 0, 0, sensor_direction_inertial(1), sensor_direction_inertial(2), sensor_direction_inertial(3), Color = "m", DisplayName="", AutoScaleFactor=1, HandleVisibility="off")
end
plot_basis(gtarg.element, "f", "--", scale = 1)
quiver3(0, 0, 0, sensor_direction_inertial(1), sensor_direction_inertial(2), sensor_direction_inertial(3), DisplayName="sensor", AutoScaleFactor=1, Color = "m")
quiver3(0, 0, 0, sun_direction(1), sun_direction(2), sun_direction(3), DisplayName="sun", AutoScaleFactor=1)
[X,Y,Z]=cylinder([0 deg2rad(cone_angle)], 50);
%axis([0 1,-1 1,-.5 .5])
M=[sun_dcm, zeros([3, 1]); zeros([1, 3]), 1];
h=surf(X,Y,Z,'Parent',hgtransform('Matrix',M),'LineStyle','none','FaceAlpha',0.4);
view([30,35])
grid on
light
axis equal
xlim([-1, 1])
ylim([-1, 1])
zlim([-1, 1])