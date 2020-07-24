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

data.weightMass = 1;
data.weightCoM = 1;
data.weightI = 1;
                                
data = InertiaIndex(data);

% all_segments = [data.pelvis, data.thorax, data.head, ...
%                 data.right_arm, data.right_forearm, data.right_hand, ...
%                 data.left_arm, data.left_forearm, data.left_hand, ...
%                 data.right_thigh, data.right_leg, data.right_foot, ...
%                 data.left_thigh, data.left_leg, data.left_foot];
all_segments = [data.left_arm, data.left_forearm];
data.segments = all_segments;

data.nSegment = length(data.segments); data.nCardinalCoor = 3;
data.massBound = ones(data.nSegment,1) * 1; % kg
data.CoMBound = ones(data.nSegment,1) * 0.1; % m
data.inertiaBound = ones(data.nSegment,1) * 0.2;

data.gravityRotationBound = pi/16;
data.kalman_optimised_filename = ['../gravity/Solutions/Do_822_F' ...
                                    num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
                                    '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
                                    '_weightQV' num2str(1) '-' num2str(0.01) ...
                                    '_gravityRotationBound=' num2str(data.gravityRotationBound) ...
                                    '_IPOPTMA57_Q.mat'];

disp('Generating Model')
[model, data] = GenerateModel(data);
disp('Loading Kalman Filter')
[model, data] = GenerateKalmanFilter(model,data);
disp('Loading Real Data')
[model, data] = GenerateRealData(model,data);
disp('Initialize Estimation')
data = saveInitialValues(model, data);
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

N_mass = data.nSegment;
w0 = [w0; data.initialMass];

N_CoM = data.nSegment * data.nCardinalCoor;
w0 = [w0; reshape(data.initialCoM',[N_CoM,1])];

N_I = data.nSegment * data.nCardinalCoor;
w0 = [w0; reshape(data.initialInertia',[N_I,1])];

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

q_opt = nan(model.nq,data.Nint+1);
v_opt = nan(model.nq,data.Nint+1);
u_opt = nan(model.nu,data.Nint);
w_opt = full(sol.x);

N_extras = N_mass + N_CoM + N_I;
for i=1:model.nq
    q_opt(i,:) = w_opt(i:model.nx+model.nu:end - N_extras)';
    v_opt(i,:) = w_opt(i+model.nq:model.nx+model.nu:end - N_extras)';
end
for i=1:model.nu
    u_opt(i,:) = w_opt(i+model.nx:model.nx+model.nu:end - N_extras)';
end
mass_opt = w_opt(end - N_extras + 1:end - N_extras + N_mass);

CoM_opt = reshape(w_opt(end - N_extras + N_mass + 1:end - N_I),[data.nCardinalCoor, data.nSegment])';

I_opt = reshape(w_opt(end - N_I + 1:end),[data.nCardinalCoor, data.nSegment])';

data.mass_opt = mass_opt;
data.CoM_opt = CoM_opt;
data.I_opt = I_opt;

data.q_opt = q_opt;
data.v_opt = v_opt;
data.u_opt = u_opt;

for i=1:length(data.segments)
    model.I{data.segments(i)} = mcI(data.mass_opt(i),data.CoM_opt(i,:), diag(data.I_opt(i,:)));
end

% disp('Calculating Simulation')
% [model, data] = GenerateSimulation(model, data);
disp('Calculating Momentum')
data = CalculateMomentum(model, data);

stats = solver.stats;
save(['Solutions/Do_822_F' num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
      '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
      '_weightQV' num2str(data.weightQV(1)) '-' num2str(data.weightQV(2)) ...
      '_bounds=' strjoin(strsplit(num2str(data.massBound'), ' ', 'CollapseDelimiters', true),',') ...
      '_' strjoin(strsplit(num2str(data.CoMBound'), ' ', 'CollapseDelimiters', true),',') ...
      '_' strjoin(strsplit(num2str(data.inertiaBound'), ' ', 'CollapseDelimiters', true),',') ...
      '_IPOPTMA57_Q.mat'],'model','data','stats')
% GeneratePlots(model, data);
% AnimatePlot(model, data, 'sol', 'kalman');

disp('mass_opt')
disp(data.mass_opt)
disp('initialMass')
disp(data.initialMass)

disp('CoM_opt')
disp(data.CoM_opt)
disp('initialCoM')
disp(data.initialCoM)

disp('I_opt')
disp(data.I_opt)
disp('initialInertia')
disp(data.initialInertia)

toc
% showmotion(model, 0:data.Duration/data.Nint:data.Duration, q_opt(:,:))
