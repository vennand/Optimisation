% Script to optimize a trajectory with 42 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
tic
run('../startup.m')
import casadi.*

data.nDoF = 42;

data.Nint = 50;% number of control nodesl
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.gravity = [0; 0; -9.81];
data.gravityRotationBound = pi/16;
data.nCardinalCoor = 3;

data.dataFile = '../data/Do_822_contact_2.c3d';
data.kalmanDataFile_q = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_Q.mat';
data.kalmanDataFile_v = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_V.mat';
data.kalmanDataFile_a = '../data/Do_822_contact_2_MOD200.00_GenderF_DoCig_A.mat';

% Spécific à Do_822_contact_2.c3d
% Le saut est entre les frames 3050 et 3386
% data.frames = 3078:3368; % Sans contact avec la trampoline
% data.frames = 3100:3311; % Sans contact avec la trampoline, interval plus sévère
data.frames = 3100:3300;
data.labels = 1:95;

data.realNint = length(data.frames);

data = adjust_number_of_interval(data);

data.weightU = 1e-7;
data.weightX = 1;
data.weightQV = [1; 0.01];

disp('Generating Model')
[model, data] = GenerateModel(data);
% disp('Loading Kalman Filter')
% [model, data] = GenerateKalmanFilter(model,data);

% Load simulated data
simulation = load('../simulation/Simulations/Do_822_simN200_simV0.001_U0.01_N50_G1,0.5,-9.7461_CoM0,0,0.17429_I0.39386,0.43611,0.18255.mat');
sim_data = simulation.data;
sim_model = simulation.model;
data.sim_data = sim_data;
data.sim_model = sim_model;

data.kalman_q = sim_data.gaussianNoiseX(1:model.nq,:);
data.kalman_v = sim_data.gaussianNoiseX(model.nq+1:end,:);
data.kalman_tau = sim_data.u;

disp('Loading Real Data')
[model, data] = GenerateRealData(model,data);
disp('Calculating Estimation')
[prob, lbw, ubw, lbg, ubg, objFunc, conFunc, objGrad, conGrad] = GenerateEstimation_Q_multiple_shooting(model, data);

% [lbw, ubw] = GenerateInitialConstraints(model, data, lbw, ubw);
% [lbw, ubw] = GenerateFinalConstraints(model, data, lbw, ubw);

options = struct;
options.ipopt.max_iter = 3000;
options.ipopt.print_level = 5;
options.ipopt.linear_solver = 'ma57';

options.ipopt.tol = 1e-6; % default: 1e-08
% options.ipopt.acceptable_tol = 1e-4; % default: 1e-06
options.ipopt.constr_viol_tol = 0.001; % default: 0.0001
% options.ipopt.acceptable_constr_viol_tol = 0.1; % default: 0.01

disp('Generating Solver')
% solver = nlpsol('solver', 'snopt', prob, options); % FAIRE MARCHER ÇA
solver = nlpsol('solver', 'ipopt', prob, options);

w0=[];
for k=1:data.Nint
    w0 = [w0; data.kalman_q(:,k); data.kalman_v(:,k)];
    w0 = [w0; data.kalman_tau(:,k)];
end
w0 = [w0; data.kalman_q(:,data.Nint+1); data.kalman_v(:,data.Nint+1)];

N_G = data.nCardinalCoor;
w0 = [w0; data.gravity];

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

q_opt = nan(model.nq,data.Nint+1);
v_opt = nan(model.nq,data.Nint+1);
u_opt = nan(model.nu,data.Nint);
w_opt = full(sol.x);

for i=1:model.nq
    q_opt(i,:) = w_opt(i:model.nx+model.nu:end - N_G)';
    v_opt(i,:) = w_opt(i+model.nq:model.nx+model.nu:end - N_G)';
end
for i=1:model.nu
    u_opt(i,:) = w_opt(i+model.nx:model.nx+model.nu:end - N_G)';
end
G_opt = w_opt(end - N_G + 1:end);

data.G_opt = G_opt;

data.q_opt = q_opt;
data.v_opt = v_opt;
data.u_opt = u_opt;

model.gravity = data.G_opt;

disp('Calculating Simulation')
[model, data] = GenerateSimulation(model, data);
disp('Calculating Momentum')
data = CalculateMomentum(model, data);

stats = solver.stats;
save(['Simulations/Do_822_F' num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
      '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
      '_weightQV' num2str(data.weightQV(1)) '-' num2str(data.weightQV(2)) ...
      '_gravityRotationBound=' num2str(data.gravityRotationBound) ...
      '_IPOPTMA57_Q.mat'],'model','data','stats')
% GeneratePlots(model, data);
% AnimatePlot(model, data, 'sol', 'kalman');
toc
% showmotion(model, 0:data.Duration/data.Nint:data.Duration, q_opt(:,:))
