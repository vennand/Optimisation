function [prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_Q_multiple_shooting(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes
dN = T/Nint;

weightQV = vertcat(data.weightQV(1) * ones(model.nq,1), data.weightQV(2) * ones(model.nq,1));

if data.optimiseInertia
    mass = zeros(4,1);
    CoM = zeros(4,3);
    I = zeros(4,3,3);
    
    pelvis = 6;
    thorax = 9;
    right_thigh = 33;
    left_thigh = 39;

    [mass(1),CoM(1,:),I(1,:,:)] = mcI(model.I{pelvis});
    [mass(2),CoM(2,:),I(2,:,:)] = mcI(model.I{thorax});
    [mass(3),CoM(3,:),I(3,:,:)] = mcI(model.I{right_thigh});
    [mass(4),CoM(4,:),I(4,:,:)] = mcI(model.I{left_thigh});
    
    data.initialMass = mass;
    data.initialCoM = CoM;
end

tau_base = SX.zeros(6,1);
x = SX.sym('x', model.nx);
u = SX.sym('u', model.nu);

L = @(x)data.weightX * ((weightQV.*x)' * (weightQV.*x));
S = @(u)data.weightU * (u'*u);

if data.optimiseGravity && data.optimiseInertia
    G = SX.sym('G',3);
    mass = SX.sym('G',4);
    CoM = SX.sym('G',4,3);
    forDyn = @(x,u,G)[  x(model.idx_v)
        FDab_Casadi_gravity_inertia( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u), G, mass, CoM, I )];
    f = Function('f', {x, u, G, mass, CoM}, {forDyn(x,u,G,mass,CoM)});
elseif data.optimiseGravity
    G = SX.sym('G',3);
    forDyn = @(x,u,G)[  x(model.idx_v)
        FDab_Casadi_gravity( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u), G )];
    f = Function('f', {x, u, G}, {forDyn(x,u,G)});
elseif data.optimiseInertia
    mass = SX.sym('G',4);
    CoM = SX.sym('G',4,3);
    forDyn = @(x,u,G)[  x(model.idx_v)
        FDab_Casadi_inertia( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u), mass, CoM, I )];
    f = Function('f', {x, u, mass, CoM}, {forDyn(x,u,mass,CoM)});
else
    forDyn = @(x,u)[  x(model.idx_v)
        FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )];
    f = Function('f', {x, u}, {forDyn(x,u)});
end

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

if data.optimiseGravity
    G = MX.sym('G',3);
end
if data.optimiseInertia
    mass = MX.sym('G',4);
    CoM = MX.sym('G',4,3);
end

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
    if data.optimiseGravity && data.optimiseInertia
        for j=1:M
            k1 = f(Xk, Uk, G, mass, CoM);
            k2 = f(Xk + DT/2 * k1, Uk, G, mass, CoM);
            k3 = f(Xk + DT/2 * k2, Uk, G, mass, CoM);
            k4 = f(Xk + DT * k3, Uk, G, mass, CoM);

            Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
        end
    elseif data.optimiseGravity
        for j=1:M
            k1 = f(Xk, Uk, G);
            k2 = f(Xk + DT/2 * k1, Uk, G);
            k3 = f(Xk + DT/2 * k2, Uk, G);
            k4 = f(Xk + DT * k3, Uk, G);

            Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
        end
    elseif data.optimiseInertia
        for j=1:M
            k1 = f(Xk, Uk, mass, CoM);
            k2 = f(Xk + DT/2 * k1, Uk, mass, CoM);
            k3 = f(Xk + DT/2 * k2, Uk, mass, CoM);
            k4 = f(Xk + DT * k3, Uk, mass, CoM);

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
    
    Jx = {Jx{:}, fJx(X_kalman(:,k+2) - Xk)};
    
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
    lbg = [lbg; -1e-12];
    ubg = [ubg;  1e-12];
end
if data.optimiseInertia
    mass_bounds = data.initialMass .* data.inertiaRelativeBound;
    
    w = {w{:}, mass};
    lbw = [lbw; data.initialMass - mass_bounds];
    ubw = [ubw; data.initialMass + mass_bounds];
    
    [N_segment, N_cardinal_coor] = size(data.initialCoM);
    CoM_bounds = zeros(size(data.initialCoM));
    for i = 1:N_cardinal_coor
        CoM_bounds(:,i) = data.initialCoM(:,i) .* data.inertiaRelativeBound;
    end
    
    for j = 1:N_segment
        w = {w{:}, CoM(j,:)'};
        lbw = [lbw; data.initialCoM(j,:)' - CoM_bounds(j,:)'];
        ubw = [ubw; data.initialCoM(j,:)' + CoM_bounds(j,:)'];
    end
end

Jx = vertcat(Jx{:});
w = vertcat(w{:});
g = vertcat(g{:});
prob = struct('f', sum(Jx)+Ju, 'x', w, 'g', g);

if nargout > 5
    objFunc = Function('J',  {w}, {Jx, Ju});
    conFunc = Function('g',  {w}, {g});
    objGrad = Function('dJ', {w}, {jacobian(Jx,w)});
    conGrad = Function('dg', {w}, {jacobian(g,w)});
end 

end