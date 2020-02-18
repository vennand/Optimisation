function [prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_multiple_shooting(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes
dN = T/Nint;

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

tau_base = SX.zeros(6,1);
x = SX.sym('x', model.nx);
u = SX.sym('u', model.nu);
markers = SX.sym('markers', N_cardinal_coor, N_markers);
is_nan  = SX.sym('is_nan', N_cardinal_coor, N_markers);

L = @(x)base_referential_coor(model, x(1:model.NB)); % Estimated marker positions, not objective function
S = @(u)data.weightU * (u'*u);

if data.optimiseGravity
    G = SX.sym('G',3);
    forDyn = @(x,u,G)[  x(model.idx_v)
        FDab_Casadi_gravity( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u), G )];
    f = Function('f', {x, u, G}, {forDyn(x,u,G)});
else
    forDyn = @(x,u)[  x(model.idx_v)
        FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )];
    f = Function('f', {x, u}, {forDyn(x,u)});
end

fJmarkers = Function('fJ', {x, markers, is_nan}, {data.weightPoints * objective_func(model,markers,is_nan,L(x))});
fJu = Function('fJu', {u}, {S(u)});

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
Jmarkers = {};
Ju = 0;
g={};
lbg = [];
ubg = [];

if data.optimiseGravity
    G = MX.sym('G',3);
end

Xk = MX.sym(['X_' '0'], model.nx);
w = {w{:}, Xk};
lbw = [lbw; model.xmin];
ubw = [ubw; model.xmax];

Jmarkers = {Jmarkers{:}, fJmarkers(Xk, markers(:,:,1), is_nan(:,:,1))};

M = 4;
DT = dN/M;
for k=0:Nint-1
    Uk = MX.sym(['U_' num2str(k)], model.nu);
    w = {w{:}, Uk};
    lbw = [lbw; model.umin];
    ubw = [ubw; model.umax];
    
%     Xk = RK4('x0',Xk,'p',Uk);
%     Xk = Xk.xf;
    if data.optimiseGravity
        for j=1:M
            k1 = f(Xk, Uk, G);
            k2 = f(Xk + DT/2 * k1, Uk, G);
            k3 = f(Xk + DT/2 * k2, Uk, G);
            k4 = f(Xk + DT * k3, Uk, G);

            Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
        end
    else
        for j=1:M
            k1 = f(Xk, Uk);
            k2 = f(Xk + DT/2 * k1, Uk);
            k3 = f(Xk + DT/2 * k2, Uk);
            k4 = f(Xk + DT * k3, Uk);

            Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
        end
    end
    
    Xkend = Xk;
    
    Ju = Ju + fJu(Uk);
    
    Xk = MX.sym(['X_' num2str(k+1)], model.nx);
    w = {w{:}, Xk};
    lbw = [lbw; model.xmin];
    ubw = [ubw; model.xmax];
    
    Jmarkers = {Jmarkers{:}, fJmarkers(Xk, markers(:,:,k+2), is_nan(:,:,k+2))};
    
    g = {g{:}, Xkend - Xk};
    lbg = [lbg; zeros(model.nx,1)];
    ubg = [ubg; zeros(model.nx,1)];
end

if data.optimiseGravity
    bounds = rotx(data.gravityRotationBound) * vertcat([0;0;0],data.gravity);
    
    w = {w{:}, G};
    lbw = [lbw; [-abs(bounds(5)); -abs(bounds(5)); data.gravity(3)]];
    ubw = [ubw; [ abs(bounds(5));  abs(bounds(5)); bounds(6)]];
    
    g = {g{:}, G'*G - data.gravity(:)' * data.gravity(:)};
    lbg = [lbg; -1e-3];
    ubg = [ubg;  1e-3];
end
if data.optimiseInertia
    bounds = data.inertiaTorsoRelativeBound;
end

Jmarkers = vertcat(Jmarkers{:});
w = vertcat(w{:});
g = vertcat(g{:});
prob = struct('f', sum(Jmarkers)+Ju, 'x', w, 'g', g);

if nargout >5
    objFunc = Function('J',  {w}, {Jmarkers, Ju});
    conFunc = Function('g',  {w}, {g});
    objGrad = Function('dJ', {w}, {jacobian(Jmarkers,w)});
    conGrad = Function('dg', {w}, {jacobian(g,w)});
end 

end
