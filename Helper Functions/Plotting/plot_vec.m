function [p] = plot_vec(vec, color, line_style, name, scale)
%PLOT_VEC Summary of this function goes here
%   Detailed explanation goes here
if name == ""
    p = quiver3([0],[0],[0],[vec(1)],[vec(2)],[vec(3)], Color = color, LineStyle = line_style, HandleVisibility = "off", AutoScaleFactor = scale);
else
    p = quiver3([0],[0],[0],[vec(1)],[vec(2)],[vec(3)], Color = color, LineStyle = line_style, DisplayName = name, AutoScaleFactor = scale);
end
end

