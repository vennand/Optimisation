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

data.gravityRotationBound = pi/16;
data.kalman_optimised_filename = ['../gravity/Solutions/Do_822_F' ...
                                    num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
                                    '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
                                    '_weightQV' num2str(1) '-' num2str(0.01) ...
                                    '_gravityRotationBound=' num2str(data.gravityRotationBound) ...
                                    '_IPOPTMA57_Q.mat'];
                                
pelvis = 6; thorax = 9; right_thigh = 33; left_thigh = 39;

data.segments = [pelvis, thorax, right_thigh, left_thigh];

data.massBound = [1; 2; 1; 1]; % kg
data.CoMBound = [0.1; 0.1; 0.1; 0.1]; % m
data.inertiaBound = [0.2; 0.2; 0.2; 0.2];
data.nSegment = 4; data.nCardinalCoor = 3;

disp('Generating Model')
[model, data] = GenerateModel(data);
disp('Loading Kalman Filter')
[model, data] = GenerateKalmanFilter(model,data);
disp('Loading Real Data')
[model, data] = GenerateRealData(model,data);
disp('Generating Simulation')
[model, data, simStateMassGrad, simStateCoMGrad, simStateInertiaGrad] = GenerateSimulation_RK4(model,data);
disp('Initialize Estimation')
data = saveInitialValues(model, data);
disp('Calculating Estimation')
[prob, lbw, ubw, lbg, ubg, ...
 objFunc, conFunc, objGrad, conGrad, ...
 stateMassGrad, stateCoMGrad, stateInertiaGrad] = ...
    GenerateEstimation_Q_multiple_shooting(model, data);
%  controlMassGrad, controlCoMGrad, controlInertiaGrad] = ...
%  JstateInertiaGrad, JcontrolInertiaGrad, JinertiaInertiaGrad] = ...

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
data.N_mass = N_mass;
w0 = [w0; data.initialMass];

N_CoM = data.nSegment * data.nCardinalCoor;
data.N_CoM = N_CoM;
w0 = [w0; reshape(data.initialCoM',[N_CoM,1])];

N_I = data.nSegment * data.nCardinalCoor;
data.N_I = N_I;
w0 = [w0; reshape(data.initialInertia',[N_I,1])];
data.w0 = w0;

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

q_opt = nan(model.nq,data.Nint+1);
v_opt = nan(model.nq,data.Nint+1);
u_opt = nan(model.nu,data.Nint);
w_opt = full(sol.x);
data.w_opt = w_opt;

N_extras = N_mass + N_CoM + N_I;
data.N_extras = N_extras;
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

disp('Calculating Simulation')
[model, data] = GenerateSimulation(model, data);
disp('Calculating Momentum')
data = CalculateMomentum(model, data);

x0 = zeros(model.nx,1);
u0 = zeros(model.nu,1);
x(model.nq+3) = 9.81/2;
x(model.nq+4) = -6;
% simState = simState(x0,u0);

data.simStateMassGrad = {};
data.simStateCoMGrad = {};
data.simStateInertiaGrad = {};
for l = 1:data.nSegment
    data.simStateMassGrad = {data.simStateMassGrad{:}, ...
            full(simStateMassGrad{l}(x0, u0, data.w0(end - N_extras + 1:end)))};
    data.simStateCoMGrad = {data.simStateCoMGrad{:}, ...
            full(simStateCoMGrad{l}(x0, u0, data.w0(end - N_extras + 1:end)))};
    data.simStateInertiaGrad = {data.simStateInertiaGrad{:}, ...
            full(simStateInertiaGrad{l}(x0, u0, data.w0(end - N_extras + 1:end)))};
end

data.objFunc_init = full(objFunc(data.w0));
data.objFunc_opt = full(objFunc(data.w_opt));
data.conFunc_init = full(conFunc(data.w0));
data.conFunc_opt = full(conFunc(data.w_opt));
data.objGrad_init = full(objGrad(data.w0));
data.objGrad_opt = full(objGrad(data.w_opt));
data.conGrad_init = full(conGrad(data.w0));
data.conGrad_opt = full(conGrad(data.w_opt));

data.stateMassGrad_init = {};
data.stateMassGrad_opt = {};
data.stateCoMGrad_init = {};
data.stateCoMGrad_opt = {};
data.stateInertiaGrad_init = {};
data.stateInertiaGrad_opt = {};

data.controlMassGrad_init = {};
data.controlMassGrad_opt = {};
data.controlCoMGrad_init = {};
data.controlCoMGrad_opt = {};
data.controlInertiaGrad_init = {};
data.controlInertiaGrad_opt = {};
for l = 1:data.nSegment
    data.stateMassGrad_init = {data.stateMassGrad_init{:}, full(stateMassGrad{l}(data.w0))};
    data.stateMassGrad_opt = {data.stateMassGrad_opt{:}, full(stateMassGrad{l}(data.w_opt))};
    data.stateCoMGrad_init = {data.stateCoMGrad_init{:}, full(stateCoMGrad{l}(data.w0))};
    data.stateCoMGrad_opt = {data.stateCoMGrad_opt{:}, full(stateCoMGrad{l}(data.w_opt))};
    data.stateInertiaGrad_init = {data.stateInertiaGrad_init{:}, full(stateInertiaGrad{l}(data.w0))};
    data.stateInertiaGrad_opt = {data.stateInertiaGrad_opt{:}, full(stateInertiaGrad{l}(data.w_opt))};
    
    data.controlMassGrad_init = {data.controlMassGrad_init{:}, full(controlMassGrad{l}(data.w0))};
    data.controlMassGrad_opt = {data.controlMassGrad_opt{:}, full(controlMassGrad{l}(data.w_opt))};
    data.controlCoMGrad_init = {data.controlCoMGrad_init{:}, full(controlCoMGrad{l}(data.w0))};
    data.controlCoMGrad_opt = {data.controlCoMGrad_opt{:}, full(controlCoMGrad{l}(data.w_opt))};
    data.controlInertiaGrad_init = {data.controlInertiaGrad_init{:}, full(controlInertiaGrad{l}(data.w0))};
    data.controlInertiaGrad_opt = {data.controlInertiaGrad_opt{:}, full(controlInertiaGrad{l}(data.w_opt))};
end

data.stateMassGrad_init_reorganized{data.nSegment} = [];
data.stateMassGrad_opt_reorganized{data.nSegment} = [];
data.stateCoMGrad_init_reorganized{data.nSegment} = [];
data.stateCoMGrad_opt_reorganized{data.nSegment} = [];
data.stateInertiaGrad_init_reorganized{data.nSegment} = [];
data.stateInertiaGrad_opt_reorganized{data.nSegment} = [];

for j = 1:data.nSegment
    for i=1:data.Nint
        data.stateMassGrad_init_reorganized{j} = [data.stateMassGrad_init_reorganized{j} ...
                                                data.stateMassGrad_init{j}(model.nx*(i-1)+1:model.nx*i)];
        data.stateMassGrad_opt_reorganized{j} = [data.stateMassGrad_opt_reorganized{j} ...
                                                data.stateMassGrad_opt{j}(model.nx*(i-1)+1:model.nx*i)];
        data.stateCoMGrad_init_reorganized{j} = [data.stateCoMGrad_init_reorganized{j} ...
                                                data.stateCoMGrad_init{j}(model.nx*(i-1)+1:model.nx*i,:)];
        data.stateCoMGrad_opt_reorganized{j} = [data.stateCoMGrad_opt_reorganized{j} ...
                                                data.stateCoMGrad_opt{j}(model.nx*(i-1)+1:model.nx*i,:)];
        data.stateInertiaGrad_init_reorganized{j} = [data.stateInertiaGrad_init_reorganized{j} ...
                                                data.stateInertiaGrad_init{j}(model.nx*(i-1)+1:model.nx*i,:)];
        data.stateInertiaGrad_opt_reorganized{j} = [data.stateInertiaGrad_opt_reorganized{j} ...
                                                data.stateInertiaGrad_opt{j}(model.nx*(i-1)+1:model.nx*i,:)];
    end
end

disp('Calculating Simulation')
[model, data, ...
    simStateMassGrad_MX, simStateCoMGrad_MX, simStateInertiaGrad_MX] = GenerateSimulation_MX(model, data);

x0 = [data.q_opt(:,1) ; data.v_opt(:,1)];

data.simStateMassGrad_MX = {};
data.simStateCoMGrad_MX = {};
data.simStateInertiaGrad_MX = {};
for l = 1:data.nSegment
    data.simStateMassGrad_MX = {data.simStateMassGrad_MX{:}, ...
            full(simStateMassGrad_MX{l}(x0, data.u_opt, data.w_opt(end - data.N_extras + 1:end)))};
    data.simStateCoMGrad_MX = {data.simStateCoMGrad_MX{:}, ...
            full(simStateCoMGrad_MX{l}(x0, data.u_opt, data.w_opt(end - data.N_extras + 1:end)))};
    data.simStateInertiaGrad_MX = {data.simStateInertiaGrad_MX{:}, ...
            full(simStateInertiaGrad_MX{l}(x0, data.u_opt, data.w_opt(end - data.N_extras + 1:end)))};
end

stats = solver.stats;
save(['Solutions/Do_822_F' num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
      '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
      '_weightQV' num2str(data.weightQV(1)) '-' num2str(data.weightQV(2)) ...
      '_massBounds=' strjoin(strsplit(num2str(data.massBound'), ' ', 'CollapseDelimiters', true),',') ...
      '_CoMBounds=' strjoin(strsplit(num2str(data.CoMBound'), ' ', 'CollapseDelimiters', true),',') ...
      '_inertiaBounds=' strjoin(strsplit(num2str(data.inertiaBound'), ' ', 'CollapseDelimiters', true),',') ...
      '_IPOPTMA57_Q.mat'],'model','data','stats')
% GeneratePlots(model, data);
% AnimatePlot(model, data, 'sol', 'kalman');
toc
% showmotion(model, 0:data.Duration/data.Nint:data.Duration, data.q_opt(:,:))
