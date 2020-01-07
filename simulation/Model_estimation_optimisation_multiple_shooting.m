% Script to optimize a trajectory with 9 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
run('startup.m')
import casadi.*

nDoF = '42';

data.Duration = 1; % Time horizon
data.Nint = 21;% number of control nodes
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.simNint = data.Nint;% number of control nodes
data.simVariance = 0.01;

data.weightU = 0.01;
data.weightPoints = 1;

disp('Generating Model')
[model, data] = GenerateModel_9(data);
disp('Generating Simulation')
[model, data] = GenerateSimulation_RK4(model,data);
disp('Calculating Estimation')
[prob, lbw, ubw, lbg, ubg] = GenerateEstimation_multiple_shooting(model, data);

[lbw, ubw] = GenerateInitialConstraints(model, data, lbw, ubw);

options = struct;
options.ipopt.max_iter = 3000;
options.ipopt.print_level = 5;

% solver = nlpsol('solver', 'snopt', prob, options); % FAIRE MARCHER ÇA
solver = nlpsol('solver', 'ipopt', prob, options);

w0=[];
for k=1:data.Nint
    w0 = [w0;  data.x(:,k)];
    w0 = [w0;  data.u(:,k)];
%     w0 = [w0;  data.x0];
%     w0 = [w0;  data.u0];
end

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

q_opt = nan(model.nq,data.Nint);
v_opt = nan(model.nq,data.Nint);
u_opt = nan(model.nu,data.Nint);
w_opt = full(sol.x);

for i=1:model.nq
    q_opt(i,:) = w_opt(i:model.nx+model.nu:end)';
    v_opt(i,:) = w_opt(i+model.nq:model.nx+model.nu:end)';
end
for i=1:model.nu
    u_opt(i,:) = w_opt(i+model.nx:model.nx+model.nu:end)';
end

GeneratePlots(model, data, q_opt, v_opt, u_opt);

% showmotion(model, 0:data.Duration/(data.Nint-1):data.Duration, q_opt(:,:))
% showmotion(model, 0:data.Duration/(data.Nint-1):data.Duration, data.xFull(1:model.nq,:))
