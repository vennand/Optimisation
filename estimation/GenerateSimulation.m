function [model, data] = GenerateSimulation(model, data)
import casadi.*

T = data.Duration; % secondes
Nint = data.Nint; % nb colloc nodes

dN = T/Nint;

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

tau_base = SX.zeros(6,1);
forDyn = @(x,u)[  x(model.idx_v)
    FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )];
x = SX.sym('x', model.nx,1);
u = SX.sym('u', model.nu,1);

f = Function('f', {x, u}, {forDyn(x,u)});

% ode = struct('x',x,'p',u,'ode',[  x(model.idx_v)
%     FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base ,u)  )]);
% opts = struct('t0',0,'tf',dN,'number_of_finite_elements',4);
% RK4 = integrator('RK4','rk',ode,opts);

Xk = [data.q_opt ; data.v_opt];
Uk = data.u_opt;
markers = zeros(N_cardinal_coor,N_markers,Nint+1);

x = Xk(:,1);
markers(:,:,1) = base_referential_coor(model, Xk(model.idx_q,1));

M = 4;
DT = dN/M;
for k=0:Nint-1
    u = Uk(:,k+1);
    
%     xRK4 = RK4('x0',x,'p',u);
%     x = xRK4.xf;
    for j=1:M
        k1 = f(x, u);
        k2 = f(x + DT/2 * k1, u);
        k3 = f(x + DT/2 * k2, u);
        k4 = f(x + DT * k3, u);

        x=x+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    
    Xk(:,k+2) = full(x);
    markers(:,:,k+2) = base_referential_coor(model, Xk(model.idx_q,k+2));
end

% Storing data
data.q_opt_sim = Xk(model.idx_q,:);
data.v_opt_sim = Xk(model.idx_v,:);
data.markers_sim = markers;

end
