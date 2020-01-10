% Script to optimize a trajectory with 42 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
run('startup.m')
import casadi.*

data.nDoF = 42;

% data.Duration = 1; % Time horizon
data.Nint = 100;% number of control nodes
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.dataFile = 'Do_822_contact_2.c3d';
data.kalmanDataFile = 'Do_822_contact_2_MOD200.00_GenderF_DoCig_Q.mat';

% Spécific à Do_822_contact_2.c3d
% Le saut est entre les frames 3050 et 3385
data.frames = 3050:3386;
data.labels = 1:95;

data.realNint = length(data.frames);

data.weightU = 0.05;
data.weightPoints = 1;

disp('Generating Model')
[model, data] = GenerateModel(data);
disp('Generating Kalman Filter')
[model, data] = GenerateKalmanFilter(model,data);
disp('Generating Real Data')
[model, data] = GenerateRealData(model,data);
disp('Calculating Estimation')
[prob, lbw, ubw, lbg, ubg, fobj] = GenerateEstimation_multiple_shooting(model, data);

% [lbw, ubw] = GenerateInitialConstraints(model, data, lbw, ubw);

options = struct;
options.ipopt.max_iter = 3000;
options.ipopt.print_level = 5;




disp('Generating Solver')
% solver = nlpsol('solver', 'snopt', prob, options); % FAIRE MARCHER ÇA
solver = nlpsol('solver', 'ipopt', prob, options);

w0=[];
for k=1:data.Nint
%     w0 = [w0;  data.x0];
    w0 = [w0;  data.kalman_q(:,k); data.x0(model.nq+1:end)];
    w0 = [w0;  data.u0];
end
w0 = [w0;  data.kalman_q(:,data.Nint+1); data.x0(model.nq+1:end)];

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

% GeneratePlots(model, data, q_opt, v_opt, u_opt);

% showmotion(model, 0:data.Duration/(data.Nint-1):data.Duration, q_opt(:,:))
