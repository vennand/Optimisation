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
markers = SX.sym('markers', N_cardinal_coor, N_markers);
is_nan  = SX.sym('is_nan', N_cardinal_coor, N_markers);

L = @(x)base_referential_coor(model, x(1:model.NB)); % Estimated marker positions, not objective function
S = @(u)data.weightU * (u'*u);

f = Function('f', {x, u}, {forDyn(x,u)});
fJu = Function('fJu', {u}, {S(u)});
fJmarkers = Function('fJ', {x, markers, is_nan}, {data.weightPoints * objective_func(model,markers,is_nan,L(x))});

% ode = struct('x',x,'p',u,'ode',[  x(model.idx_v)
%     FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )]);
% opts = struct('t0',0,'tf',dN,'number_of_finite_elements',4);
% RK4 = integrator('RK4','rk',ode,opts);

markers = data.markers;
is_nan = double(isnan(markers));

% Start with an empty NLP
w={};
lbw = [];
ubw = [];
J = 0;
g={};
lbg = [];
ubg = [];

% kalman_q = data.kalman_q;
% kalman_v = data.kalman_v;
% kalman_tau = data.kalman_tau;

Xk = MX.sym(['X_' '0'], model.nx);
% Xk = [kalman_q(:,1) ; kalman_v(:,1)];
w = {w{:}, Xk};
lbw = [lbw; model.xmin];
ubw = [ubw; model.xmax];

J = J + fJmarkers(Xk, markers(:,:,1), is_nan(:,:,1));

M = 4;
DT = dN/M;
for k=0:Nint-1
    Uk = MX.sym(['U_' num2str(k)], model.nu);
%     Uk = kalman_tau(:,k+1);
    w = {w{:}, Uk};
    lbw = [lbw; model.umin];
    ubw = [ubw; model.umax];
    
%     Xk = RK4('x0',Xk,'p',Uk);
%     Xk = Xk.xf;
    for j=1:M
        k1 = f(Xk, Uk);
        k2 = f(Xk + DT/2 * k1, Uk);
        k3 = f(Xk + DT/2 * k2, Uk);
        k4 = f(Xk + DT * k3, Uk);

        Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    Xkend = Xk;
    
    J = J + fJu(Uk);
    
    Xk = MX.sym(['X_' num2str(k+1)], model.nx);
%     Xk = [kalman_q(:,k+2) ; kalman_v(:,k+2)];
    w = {w{:}, Xk};
    lbw = [lbw; model.xmin];
    ubw = [ubw; model.xmax];
    
    J = J + fJmarkers(Xk, markers(:,:,k+2), is_nan(:,:,k+2));
    
    g = {g{:}, Xkend - Xk};
    lbg = [lbg; zeros(model.nx,1)];
    ubg = [ubg; zeros(model.nx,1)];
end

prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));

end

% Defined to be inside a CasADi function
function J = objective_func(model,markers,is_nan,estimated_markers)
J = 0;

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

n = 0;
for m = 1:N_markers
    distance_between_points = 0;
    for l = 1:N_cardinal_coor
        n = n + 1;
        distance_between_points = ...
            if_else(is_nan(n), ...
            distance_between_points, ...
            distance_between_points + (markers(n) - estimated_markers{n}).^2);
    end
    J = J + 0.5 * distance_between_points;
end
end
