% Script to optimize a trajectory with 9 DoF, 1sec time frame
% models with trapezoidal collocation
clear, clc, close all
run('startup.m')
import casadi.*

data.Duration = 1; % Time horizon
data.Nint = 20;% number of control intervals
data.odeMethod = 'rk4';%'sundials'; %'rk4';
data.obj = 'trajectory_estimation';%torque or trajectory
data.NLPMethod = 'Collocation';
data.collocMethod = 'trapezoidal';

data.simNint = 20;% number of control intervals
data.simVariance = 0.001;

[model, data] = GenerateModel('9',data);
[model, data] = GenerateSimulation_RK4(model,data);
[prob, lbw, ubw, lbg, ubg] = GenerateEstimation_colloc(model, data);

options = struct;
options.ipopt.max_iter = 3000;
options.ipopt.print_level = 5;

solver = nlpsol('solver', 'ipopt', prob, options);

w0=[];
for k=1:data.Nint
    w0 = [w0;  data.x0];
    w0 = [w0;  data.u0];
end

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw, 'lbg', lbg, 'ubg', ubg);

t_estim = linspace(data.Duration/data.Nint,data.Duration,data.Nint);
t_simu = linspace(data.Duration/data.simNint,data.Duration,data.simNint);

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

% POSITION PLOT
figure()
subplot(311)
hold on
plot(t_estim,q_opt(1:3,:),'o');
plot(t_simu,data.gaussianNoiseX(1:3,:),'x');
hold off
title('Root translation position')

subplot(312)
hold on
plot(t_estim,q_opt(4:6,:),'o');
plot(t_simu,data.gaussianNoiseX(4:6,:),'x');
hold off
title('Root rotation position')

subplot(313)
hold on
plot(t_estim,q_opt(7:9,:),'o');
plot(t_simu,data.gaussianNoiseX(7:9,:),'x');
hold off
title('Arm positions')

% VELOCITY PLOT
figure()
subplot(311)
hold on
plot(t_estim,v_opt(1:3,:),'o');
hold off
title('Root translation velocity')

subplot(312)
hold on
plot(t_estim,v_opt(4:6,:),'o');
hold off
title('Root rotation velocity')

subplot(313)
hold on
plot(t_estim,v_opt(7:9,:),'o');
hold off
title('Arm velocity')

% CONTROL PLOT
figure()
plot(t_estim,u_opt(1:3,:),'o');
title('Arms control')
% showmotion(model, data.Duration/data.Nint:data.Duration/data.Nint:data.Duration, q_opt(:,:))

% QVU = [];
% 
% for rep = 1:11
%     fprintf('***************** ITER %d **********************\n', rep)
%     if rep == 1
%         [w0, QVU, stat, feasible, ~] = optim(model, data,...
%             rep, QVU, solver, lbw, ubw, lbg, ubg, 1,...
%             [], [], struct, 'x');
%     else
%         [w0, QVU, stat, feasible, ~] = optim(model, data,...
%             rep, QVU, solver, lbw, ubw, lbg, ubg, 1,...
%             w0, feasible, stat, 'x');
%         
%     end
% end