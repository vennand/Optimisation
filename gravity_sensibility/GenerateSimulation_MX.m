function [model, data, simState, simStateGravityGrad] = GenerateSimulation_MX(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes

dN = T/Nint;

N_cardinal_coor = data.nCardinalCoor;

tau_base = SX.zeros(6,1);
x = SX.sym('x', model.nx,1);
u = SX.sym('u', model.nu,1);

G = SX.sym('G',N_cardinal_coor);
forDyn = @(x,u,G)[  x(model.idx_v)
    FDab_Casadi_gravity( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u), G )];
f = Function('f', {x, u, G}, {forDyn(x,u,G)});

% X = [data.q_opt ; data.v_opt];
x0 = MX.sym(['X_' '0'], model.nx);
u = MX.sym('U_opt', model.nu, Nint);

x = x0;

G = MX.sym('G',N_cardinal_coor);

Xk = MX.zeros(model.nx, data.Nint+1);
% markers = zeros(N_cardinal_coor,N_markers,Nint+1);

% x = X;
Xk(:,1) = x;
% markers(:,:,1) = base_referential_coor(model, Xk(model.idx_q,1));

% Xk_node = MX.zeros(size(Xk));
% markers_node = zeros(N_cardinal_coor,N_markers,Nint+1);

% x_node = x;
% Xk_node(:,1) = x_node;
% markers_node(:,:,1) = base_referential_coor(model, Xk_node(model.idx_q,1));

M = 4;
DT = dN/M;
for k=0:Nint-1
    uk = u(:,k+1);
    
    for j=1:M
        k1 = f(x, uk, G);
        k2 = f(x + DT/2 * k1, uk, G);
        k3 = f(x + DT/2 * k2, uk, G);
        k4 = f(x + DT * k3, uk, G);

        x=x+DT/6*(k1 +2*k2 +2*k3 +k4);
        
%         k1 = f(x_node, u);
%         k2 = f(x_node + DT/2 * k1, u);
%         k3 = f(x_node + DT/2 * k2, u);
%         k4 = f(x_node + DT * k3, u);
% 
%         x_node=x_node+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    Xk(:,k+2) = x;
%     markers(:,:,k+2) = base_referential_coor(model, Xk(model.idx_q,k+2));
    
%     Xk_node(:,k+2) = full(x_node);
%     markers_node(:,:,k+2) = base_referential_coor(model, Xk_node(model.idx_q,k+2));
    
%     x_node = X(:,k+2);
end

simState = Function('Xk',{x0,u,G},{Xk});

simStateGravityGrad = Function('Xk_gravity', {x0,u,G}, ...
            {jacobian(x,G)});

% Storing data
% data.q_opt_sim = Xk(model.idx_q,:);
% data.v_opt_sim = Xk(model.idx_v,:);
% data.markers_sim = markers;
% 
% data.q_opt_sim_node = Xk_node(model.idx_q,:);
% data.v_opt_sim_node = Xk_node(model.idx_v,:);
% data.markers_sim_node = markers_node;

end
