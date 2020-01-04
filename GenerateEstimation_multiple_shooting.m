function [prob, lbw, ubw, lbg, ubg] = GenerateEstimation_multiple_shooting(model, data, variables, constraints)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes

dN = T/Nint;

temp_markers = size(model.markers.coordinates);
N_cardinal_coor = temp_markers(1);
N_markers = temp_markers(2);

tau_base = SX.zeros(6,1);
forDyn = @(x,u)[  x(model.idx_v)
    FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )];
x = SX.sym('x', model.nx,1);
u = SX.sym('u', model.nu,1);
markers = SX.sym('markers', N_cardinal_coor * N_markers);

switch data.obj
    case 'twist', L = 0.5* (u'*u);
    case 'twistPond'
        if strcmpi(model.name, '10'), L = 10*(u([1 3])'*u([1 3]))+0.01*(u([2 4])'*u([2 4]));
        else, L = 0.01*(u'*u);
        end
    case 'torque', L = 0.5* (u'*u);
    case 'trajectory_estimation', L = @(x)base_referential_coor(model, x(1:model.NB)); % Estimated marker positions, not objective function
end
S = @(u)0.01* (u'*u);

f = Function('f', {x, u}, {forDyn(x,u)});
fJ = Function('fJ', {x, u, markers}, {S(u) + objective_func(model,markers,L(x))});

markers = data.gaussianNoiseMarkers;

% Start with an empty NLP
w={};
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

Xk = MX.sym(['X_' '1'], model.nx);
w = {w{:}, Xk};
lbw = [lbw; model.xmin];
ubw = [ubw; model.xmax];

Uk = MX.sym(['U_' '1'], model.nu);
w = {w{:}, Uk};
lbw = [lbw; model.umin];
ubw = [ubw; model.umax];

J = J + fJ(Xk, Uk, markers(1,:));

M = 4;
DT = dN/M;
for k=1:Nint-1
    for j=1:M
        k1 = f(Xk, Uk);
        k2 = f(Xk + DT/2 * k1, Uk);
        k3 = f(Xk + DT/2 * k2, Uk);
        k4 = f(Xk + DT * k3, Uk);

        Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    Xkend = Xk;
    
    Xk = MX.sym(['X_' num2str(k+1)], model.nx);
    w = {w{:}, Xk};
    lbw = [lbw; model.xmin];
    ubw = [ubw; model.xmax];
    
    Uk = MX.sym(['U_' num2str(k+1)], model.nu);
    w = {w{:}, Uk};
    lbw = [lbw; model.umin];
    ubw = [ubw; model.umax];
    
    g = {g{:}, Xkend - Xk};
    lbg = [lbg; zeros(model.nx,1)];
    ubg = [ubg; zeros(model.nx,1)];
    
    J = J + fJ(Xk, Uk, markers(k+1,:));
end

prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));

end


function J = objective_func(model,markers,estimated_markers)
J = 0;

temp_markers = size(model.markers.coordinates);
N_cardinal_coor = temp_markers(1);
N_markers = temp_markers(2);

n = 0;
for m = 1:N_markers
    distance_between_points = 0;
    for l = 1:N_cardinal_coor
        n = n + 1;
        distance_between_points = distance_between_points + (markers(n) - estimated_markers{n}).^2;
    end
    J = J + 0.5 * distance_between_points;
end
end
% 
% function new_value = range_conversion(old_value, old_max, old_min, new_max, new_min)
% old_range = (old_max - old_min);
% new_range = (new_max - new_min);
% new_value = (((old_value - old_min) * new_range) / old_range) + new_min;
% end
