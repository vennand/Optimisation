% Script to optimize a trajectory with 42 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
tic
run('../startup.m')
import casadi.*

data.nDoF = 42;

data.Nint = 5;% number of control nodes
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.dataFile = '../data/Do_822_contact_2.c3d';
data.kalmanDataFile_q = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_Q.mat';
data.kalmanDataFile_v = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_V.mat';
data.kalmanDataFile_a = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_A.mat';

% Spécific à Do_822_contact_2.c3d
% Le saut est entre les frames 3050 et 3386
% data.frames = 3078:3368; % Sans contact avec la trampoline
data.frames = 3220:3225;
data.labels = 1:95;

data.realNint = length(data.frames);

[data] = adjust_number_of_interval(data);

data.weightU = 10^-7;
data.weightPoints = 1;

disp('Generating Model')
[model, data] = GenerateModel(data);
disp('Loading Kalman Filter')
[model, data] = GenerateKalmanFilter(model,data);
disp('Loading Real Data')
[model, data] = GenerateRealData(model,data);
disp('Calculating Estimation')
tic
[prob, lbw, ubw, lbg, ubg] = GenerateEstimation_multiple_shooting(model, data);
toc

% [lbw, ubw] = GenerateInitialConstraints(model, data, lbw, ubw);
% [lbw, ubw] = GenerateFinalConstraints(model, data, lbw, ubw);

options = struct;
options.ipopt.max_iter = 3000;
options.ipopt.print_level = 5;
options.ipopt.linear_solver = 'ma57';

options.ipopt.tol = 10^-5; % default: 10-08
% options.ipopt.acceptable_tol = 10^-4; % default: 10-06
options.ipopt.constr_viol_tol = 0.001; % default: 0.0001
% options.ipopt.acceptable_constr_viol_tol = 0.1; % default: 0.01

% options.ipopt.hessian_approximation = 'limited-memory';

disp('Generating Solver')
% solver = nlpsol('solver', 'snopt', prob, options); % FAIRE MARCHER ÇA
solver = nlpsol('solver', 'ipopt', prob, options);

w0=[];
for k=1:data.Nint
%     w0 = [w0;  data.x0];
    w0 = [w0;  data.kalman_q(:,k); data.kalman_v(:,k)];
%     w0 = [w0;  data.u0];
    w0 = [w0;  data.kalman_tau(:,k)];
end
% w0 = [w0;  data.x0];
w0 = [w0;  data.kalman_q(:,data.Nint+1); data.kalman_v(:,data.Nint+1)];

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

stats = solver.stats;
save(['Solutions/Do_822_F' num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
      '_U' num2str(data.weightU) '_IPOPTMA57.mat'],'model','data','q_opt','v_opt','u_opt','stats')

% GeneratePlots(model, data, q_opt, v_opt, u_opt);
toc
% showmotion(model, 0:data.Duration/data.Nint:data.Duration, q_opt(:,:))