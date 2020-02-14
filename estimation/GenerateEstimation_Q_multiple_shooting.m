function [prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_Q_multiple_shooting(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes
dN = T/Nint;

if data.optimiseGravity
    G = SX.sym('G',3);
    model.gravity = G;
end

weightQV = vertcat(data.weightQV(1) * ones(model.nq,1), data.weightQV(2) * ones(model.nq,1));

tau_base = SX.zeros(6,1);
forDyn = @(x,u)[  x(model.idx_v)
    FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )];
x = SX.sym('x', model.nx);
u = SX.sym('u', model.nu);

L = @(x)data.weightX * ((weightQV.*x)' * (weightQV.*x));
S = @(u)data.weightU * (u'*u);

f = Function('f', {x, u}, {forDyn(x,u)});
fJx = Function('fJx', {x}, {L(x)});
fJu = Function('fJu', {u}, {S(u)});

% ode = struct('x',x,'p',u,'ode',[  x(model.idx_v)
%     FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )]);
% opts = struct('t0',0,'tf',dN,'number_of_finite_elements',4);
% RK4 = integrator('RK4','rk',ode,opts);

% Start with an empty NLP
w={};
lbw = [];
ubw = [];
Jx = {};
Ju = 0;
g={};
lbg = [];
ubg = [];

kalman_q = data.kalman_q;
kalman_v = data.kalman_v;
kalman_tau = data.kalman_tau;

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

if data.optimiseGravity
    bounds = rotx(data.gravityRotation) * vertcat([0;0;0],data.gravity);
    
    w = {w{:}, G};
    lbw = [lbw; [-abs(bounds(5)); -abs(bounds(5)); data.gravity(3)]];
    ubw = [ubw; [ abs(bounds(5));  abs(bounds(5)); bounds(6)]];
    
    g = {g{:}, norm(G) - norm(data.gravity)};
    lbg = [lbg; -1e-2];
    ubg = [ubg;  1e-2];
end

Jx = vertcat(Jx{:});
w = vertcat(w{:});
g = vertcat(g{:});
prob = struct('f', sum(Jx)+Ju, 'x', w, 'g', g);

if nargout >5
    objFunc = Function('J',  {w}, {Jx, Ju});
    conFunc = Function('g',  {w}, {g});
    objGrad = Function('dJ', {w}, {jacobian(Jx,w)});
    conGrad = Function('dg', {w}, {jacobian(g,w)});
end 

end
