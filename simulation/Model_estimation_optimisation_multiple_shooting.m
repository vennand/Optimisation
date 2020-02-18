% Script to optimize a trajectory with 9 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
run('startup.m')
import casadi.*

data.nDoF = 9;

data.Duration = 1; % Time horizon
data.Nint = 3;% number of control nodes
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.simNint = data.Nint;% number of control nodes
data.simVariance = 0.01;

data.weightU = 0.01;
data.weightPoints = 1;

disp('Generating Model')
[model, data] = GenerateModel(data);
disp('Generating Simulation')
[model, data] = GenerateSimulation_RK4(model,data);
disp('Calculating Estimation')
[prob, lbw, ubw, lbg, ubg] = GenerateEstimation_multiple_shooting(model, data);

[lbw, ubw] = GenerateInitialConstraints(model, data, lbw, ubw);

options = struct;
options.ipopt.max_iter = 3000;
options.ipopt.print_level = 5;

disp('Generating Solver')
% solver = nlpsol('solver', 'snopt', prob, options); % FAIRE MARCHER ÇA
solver = nlpsol('solver', 'ipopt', prob, options);

w0=[];
for k=1:data.Nint
    w0 = [w0;  data.x(:,k)];
    w0 = [w0;  data.u(:,k)];
%     w0 = [w0;  data.x0];
%     w0 = [w0;  data.u0];
end
w0 = [w0;  data.x(:,data.Nint+1)];

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

q_opt = nan(model.nq,data.Nint+1);
v_opt = nan(model.nq,data.Nint+1);
u_opt = nan(model.nu,data.Nint);
w_opt = full(sol.x);

for i=1:model.nq
    q_opt(i,:) = w_opt(i:model.nx+model.nu:end)';
    v_opt(i,:) = w_opt(i+model.nq:model.nx+model.nu:end)';
end
for i=1:model.nu
    u_opt(i,:) = w_opt(i+model.nx:model.nx+model.nu:end)';
end

data.q_opt = q_opt;
data.v_opt = v_opt;
data.u_opt = u_opt;

% GeneratePlots(model, data);
% CalculateMomentum(model, data);

% showmotion(model, 0:data.Duration/data.Nint:data.Duration, q_opt(:,:))
% showmotion(model, 0:data.Duration/data.Nint:data.Duration, data.xFull(1:model.nq,:))
