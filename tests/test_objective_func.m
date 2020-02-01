run('../startup.m')
import casadi.*

ref_rotation = refential_matrix();

tau_base = SX.zeros(6,1);
forDyn = @(x,u)[  x(model.idx_v)
    FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )];
x = SX.sym('x', model.nx,1);
u = SX.sym('u', model.nu,1);
markers = SX.sym('markers', N_cardinal_coor, N_markers);
is_nan  = SX.sym('is_nan', N_cardinal_coor, N_markers);

L = @(x)ref_rotation*base_referential_coor(model, x(1:model.NB)); % Estimated marker positions, not objective function

fJmarkers = Function('fJ', {x, markers, is_nan}, {data.weightPoints * objective_func(model,markers,is_nan,L(x))});

kalman_q = data.kalman_q;
kalman_v = data.kalman_v;

x = [kalman_q(:,1) ; kalman_v(:,1)];

estimated_markers = ref_rotation*base_referential_coor(model, x(1:model.NB));

markers = data.markers;
is_nan = double(isnan(markers));

%Reposition the referential to the model
for k=1:size(markers,3)
    markers(:,:,k) = ref_rotation * markers(:,:,k);
end

J = fJmarkers(x, markers(:,:,1), is_nan(:,:,1));
disp(J)

% % Defined to be inside a CasADi function
% function J = objective_func(model,markers,is_nan,estimated_markers)
% J = 0;
% 
% [N_cardinal_coor, N_markers] = size(model.markers.coordinates);
% 
% n = 0;
% for m = 1:N_markers
%     distance_between_points = 0;
%     for l = 1:N_cardinal_coor
%         n = n + 1;
%         distance_between_points = ...
%             if_else(is_nan(n), ...
%             distance_between_points, ...
%             distance_between_points + (markers(n) - estimated_markers{n}).^2);
%     end
%     J = J + 0.5 * distance_between_points;
% end
% end

% Defined to be inside a CasADi function
function J = objective_func(model,markers,is_nan,estimated_markers)
J = 0;

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

% n = 0;
for m = 1:N_markers
    distance_between_points = 0;
    for l = 1:N_cardinal_coor
%         n = n + 1;
        distance_between_points = ...
            if_else(is_nan(l,m), ...
            distance_between_points, ...
            distance_between_points + (markers(l,m) - estimated_markers(l,m)).^2);
    end
    J = J + 0.5 * distance_between_points;
end
end