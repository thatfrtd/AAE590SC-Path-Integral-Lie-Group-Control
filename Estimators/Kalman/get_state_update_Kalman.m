function [xtilde_k, Ptilde_k] = get_state_update_Kalman(x_k, P0, u_k, filter, options)
%GET_STATE_UPDATE_KALMAN Summary of this function goes here
%   Detailed explanation goes here
arguments (Input)
    x_k % True state
    P0 % Initial covariance 
    u_k 
    filter % Needs A_k, B_k, c_k, C_k, D_k, L_k
    options.v_k = randn(size(x_k, 1:2))
end

xtilde_k = zeros(size(x_k));
xtilde_k(:, 1) = x_k(:, 0) + chol(P0, "lower") * randn(size(x_k(:, 1))); % Sample starting estimate

Ptilde_k = zeros([size(x_k, 1), size(x_k, 1), size(x_k, 2)]);

for k = 1 : size(x_k, 2) - 1
    %% Estimated state
    xtilde_k(:, k + 1) = filter.A_k * xtilde_k(:, k) ...
                       + filter.B_k * u_k(:, k) ...
                       + filter.c_k;
    Ptilde_k(:, :, k + 1) = covariance_time_update(filter.A_k, Ptilde_k(:, :, k), filter.G_k);
    
    %% Measurement update
    % Get mesurement using true state
    y_k = filter.f_0(t_k(k + 1), x_k(:, k + 1), u_k(:, k)) + filter.g_0(x_k(:, k + 1), u_k(:, k)) * options.v_k(:, k + 1);
    
    % Get difference between expected and actual measurement
    ytilde_minus_k = innovation_process(y_k, filter.C_k, xtilde_k(:, k + 1));
    
    % Update state and covariance
    xtilde_k(:, k + 1) = estimate_measurement_update(xtilde_k(:, k + 1), filter.L_k(:, :, k + 1), ytilde_minus_k);
    Ptilde_k(:, :, k + 1) = covariance_measurement_update(filter.L_k(:, :, k + 1), filter.C_k, Ptilde_k(:, :, k + 1), filter.D_k);
end
end