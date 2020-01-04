function [prob, lbw, ubw, lbg, ubg] = GenerateEstimation_colloc(model, data, variables, constraints)
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
    case 'trajectory', L = 0.5* (qdot(7:end,1)'*qdot(7:end,1));
    case 'trajectory_estimation', L = @(x)base_referential_coor(model, x(1:model.NB)); % Estimated marker positions, not objective function
end
S = @(u)0.01* (u'*u);

f = Function('f', {x, u}, {forDyn(x,u)});
fS = Function('fJ', {x, u}, {S(u)});
fL = Function('fJ', {x, markers}, {objective_func(model,markers,L(x))});

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

fk = f(Xk,Uk);
qk = fS(Xk, Uk);
J = J + fL(Xk, markers(1,:));

for k=1:Nint-1
    Xk1 = MX.sym(['X_' num2str(k+1)], model.nx);
    w = [w, {Xk1}];
    lbw = [lbw; model.xmin];
    ubw = [ubw; model.xmax];
    
    Uk = MX.sym(['U_' num2str(k+1)], model.nu);
    w = {w{:}, Uk};
    lbw = [lbw; model.umin];
    ubw = [ubw; model.umax];
    
    fk1 = f(Xk1,Uk);
    Xkend = Xk+0.5*dN*(fk1+fk);
    
    g = {g{:}, Xkend - Xk1};
    lbg = [lbg; zeros(model.nx,1)];
    ubg = [ubg; zeros(model.nx,1)];
    
    qk1 = fS(Xk, Uk);
    J = J + 0.5*h*(qk1+qk);
    J = J + fL(Xk1, markers(k+1,:));
    
    Xk = Xk1;
    fk = fk1;
    qk = qk1;
end

prob = struct('f', J, 'x', vertcat(w{:}), 'g', vertcat(g{:}));

end

function new_value = range_conversion(old_value, old_max, old_min, new_max, new_min)
old_range = (old_max - old_min);
new_range = (new_max - new_min);
new_value = (((old_value - old_min) * new_range) / old_range) + new_min;
end

function X_e = rg4(X0,U,N,T)
M = 4; % RK4 steps per interval
DT = T/(N-1)/M;
X_e = X0;
for j=1:M
   k1 = dyn(X_e, U);
   k2 = dyn(X_e + DT/M/2 * k1, U);
   k3 = dyn(X_e + DT/M/2 * k2, U);
   k4 = dyn(X_e + DT/M * k3, U);
   X_e=X_e+DT/6*(k1 +2*k2 +2*k3 +k4);
end
end
