% Defined to be inside a CasADi function
function J = objective_func(model,markers,is_nan,estimated_markers)
import casadi.*

J = 0;

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

% n = 0;
for m = 1:N_markers
    distance_between_points = 0;
    for l = 1:N_cardinal_coor
        distance_between_points = ...
            if_else(is_nan(l,m), ...
            distance_between_points, ...
            distance_between_points + (markers(l,m) - estimated_markers(l,m)).^2);
    end
    J = J + 0.5 * distance_between_points;
end
end