% Script to optimize a trajectory with 42 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
tic
run('../startup.m')
import casadi.*

data.nDoF = 42;

data.Nint = 100;% number of control nodesl
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.simNint = 200;

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
disp('Loading Kalman Filter')
[model, data] = GenerateKalmanFilter(model,data);
disp('Loading Real Data')
[model, data] = GenerateRealData(model,data);
disp('Generating Simulation')
[model, data, simStateGravityGrad] = GenerateSimulation_RK4(model,data);
disp('Calculating Estimation')
[prob, lbw, ubw, lbg, ubg, ...
 objFunc, conFunc, objGrad, conGrad, ...
 stateGravityGrad] = ...
    GenerateEstimation_Q_multiple_shooting(model, data);

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
solver = nlpsol('solver', 'ipopt', prob, options);

w0=[];
for k=1:data.Nint
    w0 = [w0; data.kalman_q(:,k); data.kalman_v(:,k)];
    w0 = [w0; data.kalman_tau(:,k)];
end
w0 = [w0; data.kalman_q(:,data.Nint+1); data.kalman_v(:,data.Nint+1)];

N_G = data.nCardinalCoor;
data.N_G = N_G;
w0 = [w0; data.gravity];

data.w0 = w0;

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

q_opt = nan(model.nq,data.Nint+1);
v_opt = nan(model.nq,data.Nint+1);
u_opt = nan(model.nu,data.Nint);
w_opt = full(sol.x);
data.w_opt = w_opt;

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

% x0 = zeros(model.nx,1);
% u0 = zeros(model.nu,1);
% x(model.nq+3) = 9.81/2;
% x(model.nq+4) = -6;
% % simState = simState(x0,u0);
% 
% data.simStateMassGrad = {};
% data.simStateCoMGrad = {};
% data.simStateInertiaGrad = {};
% for l = 1:data.nSegment
%     data.simStateMassGrad = {data.simStateMassGrad{:}, ...
%             full(simStateMassGrad{l}(x0, u0, data.w0(end - N_extras + 1:end)))};
%     data.simStateCoMGrad = {data.simStateCoMGrad{:}, ...
%             full(simStateCoMGrad{l}(x0, u0, data.w0(end - N_extras + 1:end)))};
%     data.simStateInertiaGrad = {data.simStateInertiaGrad{:}, ...
%             full(simStateInertiaGrad{l}(x0, u0, data.w0(end - N_extras + 1:end)))};
% end

data.objFunc_init = full(objFunc(data.w0));
data.objFunc_opt = full(objFunc(data.w_opt));
data.conFunc_init = full(conFunc(data.w0));
data.conFunc_opt = full(conFunc(data.w_opt));
data.objGrad_init = full(objGrad(data.w0));
data.objGrad_opt = full(objGrad(data.w_opt));
data.conGrad_init = full(conGrad(data.w0));
data.conGrad_opt = full(conGrad(data.w_opt));

x0 = zeros(model.nx,1);
u0 = zeros(model.nu,1);
x(model.nq+3) = 9.81/2;
x(model.nq+4) = -6;

data.simStateGravityGrad = full(simStateGravityGrad(x0, u0, data.gravity));

data.stateGravityGrad_init = full(stateGravityGrad(data.w0));
data.stateGravityGrad_opt = full(stateGravityGrad(data.w_opt));

data.stateGravityGrad_init_reorganized = [];
data.stateGravityGrad_opt_reorganized = [];

for i=1:data.Nint
    data.stateGravityGrad_init_reorganized = [data.stateGravityGrad_init_reorganized ...
                                            data.stateGravityGrad_init(model.nx*(i-1)+1:model.nx*i,:)];
    data.stateGravityGrad_opt_reorganized = [data.stateGravityGrad_opt_reorganized ...
                                            data.stateGravityGrad_opt(model.nx*(i-1)+1:model.nx*i,:)];
end

disp('Calculating Simulation')
[model, data, ...
    simState, simStateGravityGrad_MX] = GenerateSimulation_MX(model, data);

x0 = [data.q_opt(:,1) ; data.v_opt(:,1)];

x_simMX = full(simState(x0, data.u_opt, data.G_opt));
data.simStateGravityGrad_MX = full(simStateGravityGrad_MX(x0, data.u_opt, data.G_opt));

stats = solver.stats;
save(['Solutions/Do_822_F' num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
      '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
      '_weightQV' num2str(data.weightQV(1)) '-' num2str(data.weightQV(2)) ...
      '_gravityRotationBound=' num2str(data.gravityRotationBound) ...
      '_IPOPTMA57_Q.mat'],'model','data','stats')
% GeneratePlots(model, data);
% AnimatePlot(model, data, 'sol', 'kalman');
toc
% showmotion(model, 0:data.Duration/data.Nint:data.Duration, data.q_opt(:,:))
