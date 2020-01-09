function [prob, lbw, ubw, lbg, ubg] = GenerateEstimation_multiple_shooting(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes

dN = T/Nint;

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

tau_base = SX.zeros(6,1);
forDyn = @(x,u)[  x(model.idx_v)
    FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )];
x = SX.sym('x', model.nx,1);
u = SX.sym('u', model.nu,1);
% markers = SX.sym('markers', N_cardinal_coor * N_markers);

L = @(x)base_referential_coor(model, x(1:model.NB)); % Estimated marker positions, not objective function
S = @(u)0.05* (u'*u);

f = Function('f', {x, u}, {forDyn(x,u)});
% fJ = Function('fJ', {x, u, markers}, {S(u) + objective_func(model,markers,L(x))});

markers = data.markers;

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

J = J + S(Uk) + objective_func(model,markers(1,:),L(Xk));
% J = J + fJ(Xk, Uk, markers(1,:));

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
    disp(['Calculating node: ', num2str(k)])
    J = J + S(Uk) + objective_func(model,markers(k+1,:),L(Xk));
%     J = J + fJ(Xk, Uk, markers(k+1,:));
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
        if class(markers(n)) == "casadi.SX"
            distance_between_points = distance_between_points + (markers(n) - estimated_markers{n}).^2;
        elseif ~isnan(markers(n)) % To deal with missing markers (maybe add check for NaN on all 3)
            distance_between_points = distance_between_points + (markers(n) - estimated_markers{n}).^2;
        end
    end
    J = J + 0.5 * distance_between_points;
end
end
