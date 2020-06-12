function [prob, lbw, ubw, lbg, ubg, ...
            objFunc, conFunc, objGrad, conGrad, ...
            stateGravityGrad] = ...
            GenerateEstimation_Q_multiple_shooting(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes
dN = T/Nint;

weightQV = vertcat(data.weightQV(1) * ones(model.nq,1), data.weightQV(2) * ones(model.nq,1));

N_cardinal_coor = data.nCardinalCoor;

tau_base = SX.zeros(6,1);
x = SX.sym('x', model.nx);
u = SX.sym('u', model.nu);

L = @(x)data.weightX * ((weightQV.*x)' * (weightQV.*x));
S = @(u)data.weightU * (u'*u);

G = SX.sym('G',N_cardinal_coor);
forDyn = @(x,u,G)[  x(model.idx_v)
    FDab_Casadi_gravity( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u), G )];
f = Function('f', {x, u, G}, {forDyn(x,u,G)});

fJx = Function('fJx', {x}, {L(x)});
fJu = Function('fJu', {u}, {S(u)});

% Start with an empty NLP
w={};
lbw = [];
ubw = [];
Jx = {};
Ju = 0;
g={};
lbg = [];
ubg = [];

G = MX.sym('G',N_cardinal_coor);

w_Xkend = {};

kalman_q = data.kalman_q;
kalman_v = data.kalman_v;

X_kalman = vertcat(kalman_q, kalman_v);

Xk = MX.sym(['X_' '0'], model.nx);
w = {w{:}, Xk};
lbw = [lbw; model.xmin];
ubw = [ubw; model.xmax];

Jx = {Jx{:}, fJx(X_kalman(:,1) - Xk)};

M = 4;
DT = dN/M;
for k=0:Nint-1
    Uk = MX.sym(['U_' num2str(k)], model.nu);
    w = {w{:}, Uk};
    lbw = [lbw; model.umin];
    ubw = [ubw; model.umax];
    
    for j=1:M
        k1 = f(Xk, Uk, G);
        k2 = f(Xk + DT/2 * k1, Uk, G);
        k3 = f(Xk + DT/2 * k2, Uk, G);
        k4 = f(Xk + DT * k3, Uk, G);

        Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    Xkend = Xk;
    
    w_Xkend = {w_Xkend{:}, Xkend};
    
    Ju = Ju + fJu(Uk);
    
    Xk = MX.sym(['X_' num2str(k+1)], model.nx);
    w = {w{:}, Xk};
    lbw = [lbw; model.xmin];
    ubw = [ubw; model.xmax];
    
    Jx = {Jx{:}, fJx(X_kalman(:,k+2) - Xk)};
    
    g = {g{:}, Xkend - Xk};
    lbg = [lbg; zeros(model.nx,1)];
    ubg = [ubg; zeros(model.nx,1)];
end

bounds = rotx(data.gravityRotationBound) * vertcat([0;0;0],data.gravity);

w = {w{:}, G};
lbw = [lbw; [-abs(bounds(5)); -abs(bounds(5)); data.gravity(3)]];
ubw = [ubw; [ abs(bounds(5));  abs(bounds(5)); bounds(6)]];

g = {g{:}, G'*G - data.gravity(:)' * data.gravity(:)};
lbg = [lbg; -1e-12];
ubg = [ubg;  1e-12];

Jx = vertcat(Jx{:});
w = vertcat(w{:});
w_Xkend = vertcat(w_Xkend{:});
g = vertcat(g{:});
prob = struct('f', sum(Jx)+Ju, 'x', w, 'g', g);

if nargout > 5
    objFunc = Function('J',  {w}, {Jx, Ju});
    conFunc = Function('g',  {w}, {g});
    objGrad = Function('dJ', {w}, {jacobian(Jx,w)});
    conGrad = Function('dg', {w}, {jacobian(g,w)});
    
    stateGravityGrad = Function('IdXk_mass', {w}, ...
            {jacobian(w_Xkend, G)});
end 

end