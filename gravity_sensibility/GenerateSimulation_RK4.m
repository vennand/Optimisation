function [model, data, simStateGravityGrad] = GenerateSimulation_RK4(model,data)
import casadi.*

simNint = data.simNint;

M = 4;
T = data.Duration; % secondes
DT = T/simNint/M; % secondes

N_cardinal_coor = data.nCardinalCoor;

tau_base = SX.zeros(6,1);
x = SX.sym('x', model.nx,1);
u = SX.sym('u', model.nu,1);

G = SX.sym('G',N_cardinal_coor);

model.gamma_q = @gamma_q;

forDyn = @(x,u,mass,CoM,I)[  x(model.idx_v)
%     FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base,u) )];
    FDgq_Casadi_gravity( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base,u), G )];

f = Function('f', {x, u, G}, {forDyn(x,u,G)});

x = MX.sym(['X_' '0'], model.nx);
u = MX.sym(['U_' '0'], model.nu);

x0 = x;
u0 = u;

G = MX.sym('G',N_cardinal_coor);

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
        k1 = f(x, u, G);
        k2 = f(x + DT/2 * k1, u, G);
        k3 = f(x + DT/2 * k2, u, G);
        k4 = f(x + DT * k3, u, G);
        
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

% data.x = Xi;
% data.u = Ui;

% simState = Function('XU',{x0,u0},{x});

simStateGravityGrad = Function('Xk_gravity', {x0,u0,G}, ...
            {jacobian(x,G)});
end