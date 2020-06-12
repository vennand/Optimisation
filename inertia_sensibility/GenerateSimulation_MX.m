function [model, data, simStateMassGrad, simStateCoMGrad, simStateInertiaGrad] = GenerateSimulation_MX(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes

dN = T/Nint;

N_cardinal_coor = data.nCardinalCoor;
N_segment = data.nSegment;

tau_base = SX.zeros(6,1);
x = SX.sym('x', model.nx,1);
u = SX.sym('u', model.nu,1);

mass = {};
CoM = {};
I = {};
for l = 1:N_segment
    mass = {mass{:}, SX.sym(['mass_' num2str(l)], 1)};
    CoM = {CoM{:}, SX.sym(['CoM_' num2str(l)], N_cardinal_coor)};
    I = {I{:}, SX.sym(['I_' num2str(l)], N_cardinal_coor)};
end
mass = vertcat(mass{:});
CoM = vertcat(CoM{:});
I = vertcat(I{:});

forDyn = @(x,u,mass,CoM,I)[  x(model.idx_v)
    FDab_Casadi_inertia( model, data, x(model.idx_q), x(model.idx_v), vertcat(tau_base,u), mass, CoM, I )];

f = Function('f', {x, u, mass, CoM, I}, {forDyn(x,u,mass,CoM,I)});

% ode = struct('x',x,'p',u,'ode',[  x(model.idx_v)
%     FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )]);
% opts = struct('t0',0,'tf',dN,'number_of_finite_elements',4);
% RK4 = integrator('RK4','rk',ode,opts);

% X = [data.q_opt ; data.v_opt];
x0 = MX.sym(['X_' '0'], model.nx);
u = MX.sym('U_opt', model.nu, Nint);

x = x0;

mass = {};
CoM = {};
I = {};
for l = 1:N_segment
    mass = {mass{:}, MX.sym(['mass_' num2str(l)], 1)};
    CoM = {CoM{:}, MX.sym(['CoM_' num2str(l)], N_cardinal_coor)};
    I = {I{:}, MX.sym(['I_' num2str(l)], N_cardinal_coor)};
end
mass = vertcat(mass{:});
CoM = vertcat(CoM{:});
I = vertcat(I{:});

% Xk = MX.zeros(model.nx, data.Nint+1);
% markers = zeros(N_cardinal_coor,N_markers,Nint+1);

% x = X;
% Xk(:,1) = x;
% markers(:,:,1) = base_referential_coor(model, Xk(model.idx_q,1));

% Xk_node = MX.zeros(size(Xk));
% markers_node = zeros(N_cardinal_coor,N_markers,Nint+1);

% x_node = X;
% Xk_node(:,1) = x_node;
% markers_node(:,:,1) = base_referential_coor(model, Xk_node(model.idx_q,1));

M = 4;
DT = dN/M;
for k=0:Nint-1
    uk = u(:,k+1);
    
%     xRK4 = RK4('x0',x,'p',u);
%     x = xRK4.xf;
    for j=1:M
        k1 = f(x, uk, mass, CoM, I);
        k2 = f(x + DT/2 * k1, uk, mass, CoM, I);
        k3 = f(x + DT/2 * k2, uk, mass, CoM, I);
        k4 = f(x + DT * k3, uk, mass, CoM, I);

        x=x+DT/6*(k1 +2*k2 +2*k3 +k4);
        
%         k1 = f(x_node, u);
%         k2 = f(x_node + DT/2 * k1, u);
%         k3 = f(x_node + DT/2 * k2, u);
%         k4 = f(x_node + DT * k3, u);
% 
%         x_node=x_node+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
%     Xk(:,k+2) = full(x);
%     markers(:,:,k+2) = base_referential_coor(model, Xk(model.idx_q,k+2));
    
%     Xk_node(:,k+2) = full(x_node);
%     markers_node(:,:,k+2) = base_referential_coor(model, Xk_node(model.idx_q,k+2));
    
%     x_node = X(:,k+2);
end

inertia_var = {};
inertia_var = {inertia_var{:}, mass};
inertia_var = {inertia_var{:}, CoM};
inertia_var = {inertia_var{:}, I};
inertia_var = vertcat(inertia_var{:});

% simState = Function('XU',{x0,u0},{x});

% simStateMassGrad = Function('dxMass', {mass}, {jacobian(x,mass)});
% simStateCoMGrad = Function('dxCoM', {CoM}, {jacobian(x,CoM)});
% simStateInertiaGrad = Function('dxI', {I}, {jacobian(x,I)});

% N_mass = N_segment;
% N_CoM = N_segment * N_cardinal_coor;
% N_I = N_segment * N_cardinal_coor;
% N_extras = N_mass + N_CoM + N_I;

simStateMassGrad = {};
simStateCoMGrad = {};
simStateInertiaGrad = {};

for l = 1:N_segment
    simStateMassGrad = {simStateMassGrad{:}, Function('IdXk_mass', {x0,u,inertia_var}, ...
        {jacobian(x,mass(l))})};
    simStateCoMGrad = {simStateCoMGrad{:}, Function('IdXk_CoM', {x0,u,inertia_var}, ...
        {jacobian(x,CoM(3*l-2:3*l))})};
    simStateInertiaGrad = {simStateInertiaGrad{:}, Function('IdXk_I', {x0,u,inertia_var}, ...
        {jacobian(x,I(3*l-2:3*l))})};
end

% Storing data
% data.q_opt_sim = Xk(model.idx_q,:);
% data.v_opt_sim = Xk(model.idx_v,:);
% data.markers_sim = markers;
% 
% data.q_opt_sim_node = Xk_node(model.idx_q,:);
% data.v_opt_sim_node = Xk_node(model.idx_v,:);
% data.markers_sim_node = markers_node;

end
