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
t_k = 0:control_delta_t:10;
u_k = 0 * ones([G.dim, numel(t_k) - 1]);

R = eye(G.dim);
S_g = 1.6e9 * eye(G.dim);
S_twist = 8e8 * eye(G.dim);

%% Define Initial Condition and Target
g0 = SO3_RotationMatrix(angle2dcm(0, deg2rad(0), deg2rad(0))); % Initial DCM
twist0 = [0; 0; 0]; % Initial angular velocity

gtarg = SO3_RotationMatrix(angle2dcm(0, deg2rad(30), deg2rad(30))'); % Target DCM
twisttarg = zeros([3, 1]);

%%
% Define initial state covariance
sigma_xhat0 = [0.1e-1; ... % _x
               1e-1; ... % _y
               0.4e-1]; ... % _z
P0 = diag(sigma_xhat0 .^ 2);

P_k = zeros([3, 3, numel(t_k)]);
P_k(:, :, 1) = P0;

filter.A_k = G.Ad(G.Exp(twist0 * control_delta_t));
filter.G_k = B(0, g0) * sigma / sqrt(noise_delta_t) * control_delta_t;
for k = 1 : numel(t_k) - 1
    P_k(:, :, k + 1) = filter.A_k * P_k(:, :, k) * filter.A_k' ...
                          + filter.G_k * filter.G_k';
end

P_fail = 0.01;
backoff_coef = sqrt(chi2inv(1 - P_fail, G.dim));

%% Sample
m = 1000;

g_k = createArray(1, m, class(g0)); % array of group elements
g_k_end = createArray(1, m, class(g0)); % array of group elements

r_array = zeros([2, m]);
for k = 1 : m
    g_k(k) = sample_right_gaussian(gtarg, P_k(:, :, 1));
    g_k_end(k) = sample_right_gaussian(gtarg, P_k(:, :, end));
end

%% Constraint
cone_angle = 50; % [deg]
sun_dcm = angle2dcm(1.4, -0.8, 2.2);
sun_direction = sun_dcm * [0; 0; 1];
sensor_direction_body = angle2dcm(0, 0, 0) * [1; 0; 0];
sensor_direction_inertial = gtarg.element * sensor_direction_body;

crossvec = -cross(sun_direction, sensor_direction_inertial);
crossvec = crossvec / norm(crossvec);

L = backoff_coef * chol(P_k(:, :, end));

crossvec = L * crossvec;

sun_constraint = @(lambda) cosd(cone_angle) - sun_direction' * gtarg.rplus(lambda).element * sensor_direction_body;
sun_constraint_approx = @(lambda) cosd(cone_angle) - sun_direction' * gtarg.element * (eye(3) + skew(lambda)) * sensor_direction_body;

%%
syms lambda [3, 1]
syms sun_direction_sym [3, 1]
syms sensor_direction_body_sym [3, 1]
syms R_sym [3, 3]
syms L_sym [3, 3]
sun_constraint_approx_sym = -sun_direction_sym.' * R_sym * (skew(L_sym * lambda)) * sensor_direction_body_sym;

%%
alpha = L' * gtarg.element' * [-(sun_direction(2) * sensor_direction_inertial(3) - sun_direction(3) * sensor_direction_inertial(2)); ...
                                 sun_direction(1) * sensor_direction_inertial(3) - sun_direction(3) * sensor_direction_inertial(1); ...
                               -(sun_direction(1) * sensor_direction_inertial(2) - sun_direction(2) * sensor_direction_inertial(1))];
alpha = L * alpha / norm(alpha);
%%

tau = randn([G.dim, 10000]);
tau = L * tau ./ vecnorm(tau);

for i = 1 : size(tau, 2)
    const(i) = sun_constraint(tau(:, i));
end

figure
scatter3(tau(1, :), tau(2, :), tau(3, :), 1, const); hold on
quiver3(0, 0, 0, alpha(1), alpha(2), alpha(3), 1.5)
%quiver3(0, 0, 0, crossvec(1), crossvec(2), crossvec(3), 1.5)
axis equal
colorbar


%%
[c, i_min] = min(const);
tau_worst = tau(:, i_min);

sun_constraint_approx(alpha)
sun_constraint(tau_worst)
sun_constraint(alpha)
sun_constraint(crossvec)


%%
figure
plot_basis(gtarg.element, "mean", "--", scale = 1);
for i = 1 : m
    %plot_basis(g_k(i).element, "", "-", scale = 0.7);
    perturb_sensor = g_k(i).element * sensor_direction_body;
    quiver3(0, 0, 0, perturb_sensor(1), perturb_sensor(2), perturb_sensor(3), AutoScaleFactor=0.9, HandleVisibility="off", Color="m", LineStyle="-")

    perturb_sensor = g_k_end(i).element * sensor_direction_body;
    quiver3(0, 0, 0, perturb_sensor(1), perturb_sensor(2), perturb_sensor(3), AutoScaleFactor=0.7, HandleVisibility="off", Color="b", LineStyle="-")
end
%perturb_sensor_worst = gtarg.element * (eye(3) + skew(alpha)) * sensor_direction_body;
perturb_sensor_worst = gtarg.rplus(alpha).element * sensor_direction_body;
quiver3(0, 0, 0, perturb_sensor_worst(1), perturb_sensor_worst(2), perturb_sensor_worst(3), AutoScaleFactor=0.7, HandleVisibility="off", Color="k", LineStyle="--")

perturb_sensor_worstest = gtarg.rplus(tau_worst).element * sensor_direction_body;
quiver3(0, 0, 0, perturb_sensor_worstest(1), perturb_sensor_worstest(2), perturb_sensor_worstest(3), AutoScaleFactor=0.7, HandleVisibility="off", Color="k", LineStyle="-")

quiver3(0, 0, 0, sun_direction(1), sun_direction(2), sun_direction(3), DisplayName="sun", AutoScaleFactor=1)
quiver3(0, 0, 0, sensor_direction_inertial(1), sensor_direction_inertial(2), sensor_direction_inertial(3), DisplayName="sensor", AutoScaleFactor=1)
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