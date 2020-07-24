function [prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_Q_multiple_shooting(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes
dN = T/Nint;

weightQV = vertcat(data.weightQV(1) * ones(model.nq,1), data.weightQV(2) * ones(model.nq,1));

N_cardinal_coor = data.nCardinalCoor;
N_segment = data.nSegment;

tau_base = SX.zeros(6,1);
x = SX.sym('x', model.nx);
u = SX.sym('u', model.nu);

L = @(x)data.weightX * ((weightQV.*x)' * (weightQV.*x));
S = @(u)data.weightU * (u'*u);
Z = @(mass, CoM, I)data.weightMass * mass'*mass + ...
                   data.weightCoM * sum(sum(CoM.*CoM)) + ...
                   data.weightI * sum(sum(I.*I));

mass = SX.sym('M',N_segment);
CoM = SX.sym('CoM',N_segment,N_cardinal_coor);
I = SX.sym('I',N_segment,N_cardinal_coor);
forDyn = @(x,u,mass,CoM,I)[  x(model.idx_v)
    FDab_Casadi_inertia( model, data, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u), mass, CoM, I )];
f = Function('f', {x, u, mass, CoM, I}, {forDyn(x,u,mass,CoM,I)});

fJx = Function('fJx', {x}, {L(x)});
fJu = Function('fJu', {u}, {S(u)});
fJI = Function('fJI', {mass, CoM, I}, {Z(mass, CoM, I)});

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

mass = MX.sym('M',N_segment);
% Étrange, mais il faut le faire comme ça, sinon l'expression finale
% n'est pas purement symbolique
% CoM = MX.sym('CoM',N_segment,N_cardinal_coor);
% L'expression précédente, quoique plus simple, ne fonctionne pas
CoM = MX.sym('CoM',N_segment*N_cardinal_coor,1);
CoM = reshape(CoM,N_cardinal_coor,N_segment)';
I = MX.sym('CoM',N_segment*N_cardinal_coor,1);
I = reshape(I,N_cardinal_coor,N_segment)';

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
    
%     Xk = RK4('x0',Xk,'p',Uk);
%     Xk = Xk.xf;
    
    for j=1:M
        k1 = f(Xk, Uk, mass, CoM, I);
        k2 = f(Xk + DT/2 * k1, Uk, mass, CoM, I);
        k3 = f(Xk + DT/2 * k2, Uk, mass, CoM, I);
        k4 = f(Xk + DT * k3, Uk, mass, CoM, I);

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

w = {w{:}, mass};
lbw = [lbw; data.initialMass - data.massBound];
ubw = [ubw; data.initialMass + data.massBound];

g = {g{:}, mass'*mass - data.initialMass(:)'*data.initialMass(:)};
lbg = [lbg; -0.1];
ubg = [ubg;  0.1];

g = {g{:}, mass};
lbg = [lbg; zeros(N_segment,1)];
ubg = [ubg;  10*ones(N_segment,1)];

CoM_lower_bounds = zeros(N_segment,N_cardinal_coor);
CoM_upper_bounds = zeros(N_segment,N_cardinal_coor);
for i = 1:N_cardinal_coor
    CoM_lower_bounds(:,i) = data.initialCoM(:,i) - data.CoMBound;
    CoM_upper_bounds(:,i) = data.initialCoM(:,i) + data.CoMBound;
end

w = {w{:}, reshape(CoM',N_segment*N_cardinal_coor,1)};
lbw = [lbw; reshape(CoM_lower_bounds',N_segment*N_cardinal_coor,1)];
ubw = [ubw; reshape(CoM_upper_bounds',N_segment*N_cardinal_coor,1)];

% g = {g{:}, reshape(CoM(:,1:2)',N_segment*(N_cardinal_coor-1),1)};
% lbg = [lbg; zeros(N_segment*(N_cardinal_coor-1),1)];
% ubg = [ubg;  zeros(N_segment*(N_cardinal_coor-1),1)];

I_lower_bounds = zeros(N_segment,N_cardinal_coor);
I_upper_bounds = zeros(N_segment,N_cardinal_coor);
for i = 1:N_cardinal_coor
    I_lower_bounds(:,i) = data.initialInertia(:,i) - data.inertiaBound;
    I_upper_bounds(:,i) = data.initialInertia(:,i) + data.inertiaBound;
end

w = {w{:}, reshape(I',N_segment*N_cardinal_coor,1)};
lbw = [lbw; reshape(I_lower_bounds',N_segment*N_cardinal_coor,1)];
ubw = [ubw; reshape(I_upper_bounds',N_segment*N_cardinal_coor,1)];

g = {g{:}, reshape(I',N_segment*N_cardinal_coor,1)};
lbg = [lbg; zeros(N_segment*N_cardinal_coor,1)];
ubg = [ubg;  ones(N_segment*N_cardinal_coor,1)];

JI = fJI(mass - data.initialMass, ...
         CoM - data.initialCoM, ...
         I - data.initialInertia);

Jx = vertcat(Jx{:});
w = vertcat(w{:});
g = vertcat(g{:});
% prob = struct('f', sum(Jx)+Ju+JI, 'x', w, 'g', g);
prob = struct('f', sum(Jx)+Ju, 'x', w, 'g', g);
% prob = struct('f', sum(Jx), 'x', w, 'g', g);

if nargout > 5
    objFunc = Function('J',  {w}, {Jx, Ju, JI});
    conFunc = Function('g',  {w}, {g});
    objGrad = Function('dJ', {w}, {jacobian(Jx,w)});
    conGrad = Function('dg', {w}, {jacobian(g,w)});
end 

end