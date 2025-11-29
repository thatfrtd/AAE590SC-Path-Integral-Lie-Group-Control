function [p] = plot_basis(n_i, basis_name, line_style, options)
arguments
    n_i
    basis_name
    line_style
    options.scale = 0.9
end
%PLOT_BASIS Summary of this function goes here
%   Detailed explanation goes here

% Plot the basis
p = [];
p.p1 = plot_vec(n_i(:, 1), "b", line_style, "$\hat " + basis_name + "_1$", options.scale); hold on;
p.p2 = plot_vec(n_i(:, 2), "r", line_style, "$\hat " + basis_name + "_2$", options.scale); hold on;
p.p3 = plot_vec(n_i(:, 3), "g", line_style, "$\hat " + basis_name + "_3$", options.scale);
legend(interpreter = "latex")
xlim([-1,1])
ylim([-1,1])
zlim([-1,1])
xlabel("X")
ylabel("Y")
zlabel("Z")

end

