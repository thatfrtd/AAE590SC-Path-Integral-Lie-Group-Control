function [] = plot_2DoF_MC_trajectories(x_mean, x_MC, t_k, P_k, options)
arguments
    x_mean
    x_MC
    t_k
    P_k
    options.x_ref_solution = []
    options.ref_sol_name = "Deterministic"
    options.title = "Monte Carlo Simulation of 2DoF Rocket Landing with and without Optimized Feedback Policy"
end
%PLOT_3DOF_MC_TRAJECTORIES Summary of this function goes here
%   Detailed explanation goes here

figure
tiledlayout(1, 1, "TileSpacing","compact")

%% Feedback controlled
nexttile

%proj_P_r = project_ellipsoid(P_k, [1,2]);

%[P_eigvecs, P_eigvals] = pageeig(proj_P_r);
[P_eigvecs, P_eigvals] = pageeig(P_k(1:2, 1:2, :)); % Looks more correct then projecting ellipsoid...

%X_k = zeros(size(P_k));
thetas = reshape(linspace(0, 2 * pi, 100), 1, []);
ellipse_3sigma = zeros([2, 100, numel(t_k)]);
for k = 1:numel(t_k)
    ellipse_3sigma(:, :, k) = x_mean(1:2, k) + P_eigvecs(:, :, k) * [3 * sqrt(P_eigvals(1, 1, k)) * cos(thetas); 3 * sqrt(P_eigvals(2, 2, k)) * sin(thetas)];
    %X_k(:, :, k) = chol(P_k(:, :, k), "lower");
end

plot(squeeze(x_MC(1, :, :)), squeeze(x_MC(2, :, :)), Color = [192, 192, 192] / 256, HandleVisibility='off'); hold on
plot(x_mean(1, :), x_mean(2, :), Color = [30, 144, 255] / 256, LineWidth=1, DisplayName="Nominal"); hold on
if ~isempty(options.x_ref_solution)
    plot(options.x_ref_solution(1, :), options.x_ref_solution(2, :), Color = "r", DisplayName = options.ref_sol_name); hold on
end
plot(squeeze(ellipse_3sigma(1, :, 2:end)), squeeze(ellipse_3sigma(2, :, 2:end)), Color = "k", HandleVisibility='off'); hold on
title("")
xlabel("X [km]")
ylabel("Y [km]")
legend(location = "best")
axis equal
plot(squeeze(ellipse_3sigma(1, :, 1)), squeeze(ellipse_3sigma(2, :, 1)), Color = "k", DisplayName="Covariance"); hold off
grid on
end

