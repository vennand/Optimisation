% Script to optimize a trajectory with 9 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
run('../startup.m')
import casadi.*

data.nDoF = 42;

data.Duration = 1; % Time horizon
data.Nint = 50;% number of control nodes
data.odeMethod = 'rk4';
data.NLPMethod = 'MultipleShooting';

data.simNint = 200;% number of control nodes
data.simVariance = 0.001;

data.weightU = 0.01;
data.weightPoints = 1;

disp('Generating Model')
[model, data] = GenerateModel(data);

model.gravity = [1; 0.5; -sqrt(9.81^2 - 1^2 -0.5^2)];
data.angle_gravity = acos(dot(model.gravity,[0;0;-9.81])/(norm(model.gravity)*norm([0;0;-9.81])))*180/pi;

data = InertiaIndex(data);

% all_segments = [data.pelvis, data.thorax, data.head, ...
%                 data.right_arm, data.right_forearm, data.right_hand, ...
%                 data.left_arm, data.left_forearm, data.left_hand, ...
%                 data.right_thigh, data.right_leg, data.right_foot, ...
%                 data.left_thigh, data.left_leg, data.left_foot];
all_segments = [data.left_arm, data.left_forearm];
data.sim_segments = all_segments;

for i = 1:length(all_segments)
    [mass_init, CoM_init, I_init] = mcI(model.I{all_segments(i)});
    rng(0,'twister');
    mass_sim = mass_init * (rand(1)/10+1);
    CoM_sim = CoM_init .* (rand(3,1)/10+1);
    I_sim = I_init .* (rand(3,3)/10+1);
    model.I{all_segments(i)} = mcI(mass_sim, CoM_sim, I_sim); %torse
end

disp('Generating Simulation')
[model, data] = GenerateSimulation_RK4(model,data);

save(['Simulations/Do_822_simN' num2str(data.simNint) ...
      '_simV' num2str(data.simVariance) ...
      '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
      '_G_' strjoin(strsplit(num2str(model.gravity'), ' ', 'CollapseDelimiters', true),',') ...
      '_Segments_' strjoin(strsplit(num2str(all_segments), ' ', 'CollapseDelimiters', true),',') ...
      '.mat'],'model','data')

q = data.xFull(1:model.nq,:);
% showmotion(model, 0:data.Duration/data.simNint:data.Duration, q(:,:))
  
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

% GeneratePlots(model, data, q_opt, v_opt, u_opt);
% CalculateMomentum(model, data);

% showmotion(model, 0:data.Duration/data.Nint:data.Duration, q_opt(:,:))
% showmotion(model, 0:data.Duration/data.Nint:data.Duration, data.xFull(1:model.nq,:))
