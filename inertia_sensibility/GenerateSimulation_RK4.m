function [model, data, simStateMassGrad, simStateCoMGrad, simStateInertiaGrad] = GenerateSimulation_RK4(model,data)
import casadi.*

simNint = data.simNint;

M = 4;
T = data.Duration; % secondes
DT = T/simNint/M; % secondes

N = model.NB;

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

model.gamma_q = @gamma_q;

forDyn = @(x,u,mass,CoM,I)[  x(model.idx_v)
%     FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base,u) )];
    FDgq_Casadi_inertia( model, data, x(model.idx_q), x(model.idx_v), vertcat(tau_base,u), mass, CoM, I )];

f = Function('f', {x, u, mass, CoM, I}, {forDyn(x,u,mass,CoM,I)});

x = MX.sym(['X_' '0'], model.nx);
u = MX.sym(['U_' '0'], model.nu);

x0 = x;
u0 = u;

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

% % Initial velocities
% x(N+1) = -model.gravity(1)/2;
% x(N+2) = -model.gravity(2)/2;
% x(N+3) = -model.gravity(3)/2;
% x(N+3) = 9.81/2;
% x(N+4) = -6;
% x(N+7) = -1.5;
% x(N+8) = 2;
% x(N+9) = 2;

% u(1) = -2;

Xi = MX.zeros(model.nx,simNint+1);
Ui = MX.zeros(model.nu,simNint);
Xi(:,1) = x;
for i = 2:simNint+1
    Ui(:,i-1) = u;
    
    for j=1:M
        k1 = f(x, u, mass, CoM, I);
        k2 = f(x + DT/2 * k1, u, mass, CoM, I);
        k3 = f(x + DT/2 * k2, u, mass, CoM, I);
        k4 = f(x + DT * k3, u, mass, CoM, I);
        
        x=x+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    Xi(:,i) = x;
    
%     if i == floor(simNint/2)
%         u(1) = 0.2;
%     elseif i == floor(simNint*3/8-1)
%         u(1) = 0.2;
%     elseif i == floor(simNint*4/8-1)
%         u(4) = 0.2;
%     elseif i == floor(simNint*5/8-1)
%         u(1) = -0.2;
%     elseif i == floor(simNint*6/8-1)
%         u(4) = -0.2;
%     elseif i == floor(simNint*7/8-1)
%         u(1) = 0.2;
%     end
    
end

% N_mass = data.nSegment;
% N_CoM = data.nSegment * data.nCardinalCoor;
% N_I = data.nSegment * data.nCardinalCoor;
% N_extras = N_mass + N_CoM + N_I;

inertia_var = {};
inertia_var = {inertia_var{:}, mass};
inertia_var = {inertia_var{:}, CoM};
inertia_var = {inertia_var{:}, I};
inertia_var = vertcat(inertia_var{:});

% CoM = reshape(CoM',N_segment*N_cardinal_coor,1);
% I = reshape(I',N_segment*N_cardinal_coor,1);

% data.x = Xi;
% data.u = Ui;

% simState = Function('XU',{x0,u0},{x});

% simStateMassGrad = Function('dxMass', {mass}, {jacobian(x,mass)});
% simStateCoMGrad = Function('dxCoM', {CoM}, {jacobian(x,CoM)});
% simStateInertiaGrad = Function('dxI', {I}, {jacobian(x,I)});

N_mass = N_segment;
N_CoM = N_segment * N_cardinal_coor;
N_I = N_segment * N_cardinal_coor;
N_extras = N_mass + N_CoM + N_I;

simStateMassGrad = {};
simStateCoMGrad = {};
simStateInertiaGrad = {};

for l = 1:N_segment
    simStateMassGrad = {simStateMassGrad{:}, Function('IdXk_mass', {x0,u0,inertia_var}, ...
        {jacobian(x,inertia_var(end - N_extras + l))})};
    simStateCoMGrad = {simStateCoMGrad{:}, Function('IdXk_CoM', {x0,u0,inertia_var}, ...
        {jacobian(x,inertia_var(end - N_extras + N_mass + 3*l-2:end - N_extras + N_mass + 3*l))})};
    simStateInertiaGrad = {simStateInertiaGrad{:}, Function('IdXk_I', {x0,u0,inertia_var}, ...
        {jacobian(x,inertia_var(end - N_I + 3*l-2:end - N_I + 3*l))})};
end
end