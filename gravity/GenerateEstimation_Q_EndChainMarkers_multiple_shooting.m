function [prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_Q_EndChainMarkers_multiple_shooting(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes
dN = T/Nint;

weightQV = vertcat(data.weightQV(1) * ones(model.nq,1), data.weightQV(2) * ones(model.nq,1));

marker_num = [37:39 59:61 76:78 93:95];

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

tau_base = SX.zeros(6,1);
x = SX.sym('x', model.nx);
u = SX.sym('u', model.nu);
markers = SX.sym('markers', N_cardinal_coor, N_markers);
is_nan  = SX.sym('is_nan', N_cardinal_coor, N_markers);

L = @(x)data.weightX * ((weightQV.*x)' * (weightQV.*x));
L_markers = @(x)base_referential_coor(model, x(1:model.NB)); % Estimated marker positions, not objective function
S = @(u)data.weightU * (u'*u);

G = SX.sym('G',N_cardinal_coor);
forDyn = @(x,u,G)[  x(model.idx_v)
    FDab_Casadi_gravity( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u), G )];
f = Function('f', {x, u, G}, {forDyn(x,u,G)});

fJx = Function('fJx', {x}, {L(x)});
fJmarkers = Function('fJ', {x, markers, is_nan}, {data.weightPoints * objective_func(model,markers,is_nan,L_markers(x))});
fJu = Function('fJu', {u}, {S(u)});

% ode = struct('x',x,'p',u,'ode',[  x(model.idx_v)
%     FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )]);
% opts = struct('t0',0,'tf',dN,'number_of_finite_elements',4);
% RK4 = integrator('RK4','rk',ode,opts);

markers = data.markers;
for n = 1:N_markers
    if ~ismember(n,marker_num)
        markers(:,n,:) = NaN;
    end
end
is_nan = double(isnan(markers));

% Start with an empty NLP
w={};
lbw = [];
ubw = [];
Jx = {};
Jmarkers = {};
Ju = 0;
g={};
lbg = [];
ubg = [];

G = MX.sym('G',N_cardinal_coor);

kalman_q = data.kalman_q;
kalman_v = data.kalman_v;

X_kalman = vertcat(kalman_q, kalman_v);

Xk = MX.sym(['X_' '0'], model.nx);
w = {w{:}, Xk};
lbw = [lbw; model.xmin];
ubw = [ubw; model.xmax];

Jx = {Jx{:}, fJx(X_kalman(:,1) - Xk)};
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
    for j=1:M
        k1 = f(Xk, Uk, G);
        k2 = f(Xk + DT/2 * k1, Uk, G);
        k3 = f(Xk + DT/2 * k2, Uk, G);
        k4 = f(Xk + DT * k3, Uk, G);

        Xk=Xk+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    Xkend = Xk;
    
    Ju = Ju + fJu(Uk);
    
    Xk = MX.sym(['X_' num2str(k+1)], model.nx);
    w = {w{:}, Xk};
    lbw = [lbw; model.xmin];
    ubw = [ubw; model.xmax];
    
    Jx = {Jx{:}, fJx(X_kalman(:,k+2) - Xk)};
    Jmarkers = {Jmarkers{:}, fJmarkers(Xk, markers(:,:,k+2), is_nan(:,:,k+2))};
    
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
Jmarkers = vertcat(Jmarkers{:});
w = vertcat(w{:});
g = vertcat(g{:});
prob = struct('f', sum(Jx)+sum(Jmarkers)+Ju, 'x', w, 'g', g);

if nargout > 5
    objFunc = Function('J',  {w}, {Jx, Jmarkers, Ju});
    conFunc = Function('g',  {w}, {g});
    objGrad = Function('dJ', {w}, {jacobian(Jx,w)});
    conGrad = Function('dg', {w}, {jacobian(g,w)});
end 

end