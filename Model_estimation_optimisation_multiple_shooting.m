% Script to optimize a trajectory with 9 or 12 DoF, 1sec time frame
% models with multiple shooting
clear, clc, close all
run('startup.m')
import casadi.*

nDoF = '9';

data.Duration = 1; % Time horizon
data.Nint = 21;% number of control nodes
data.odeMethod = 'rk4';
data.obj = 'trajectory_estimation';
data.NLPMethod = 'MultipleShooting';

data.simNint = data.Nint;% number of control nodes
data.simVariance = 0.01;

[model, data] = GenerateModel(nDoF,data);
[model, data] = GenerateSimulation_RK4(model,data);
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

t_estim = linspace(0,data.Duration,data.Nint);
t_simu = linspace(0,data.Duration,data.simNint);

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
subplot(211)
hold on
plot(t_estim,q_opt(1:3,:),'o');
plot(t_simu,data.xFull(1:3,:),'x');
plot(t_simu,data.gaussianNoiseXFull(1:3,:),'+');
hold off
title('Root translation position')

subplot(212)
hold on
plot(t_estim,q_opt(4:6,:),'o');
plot(t_simu,data.xFull(4:6,:),'x');
plot(t_simu,data.gaussianNoiseXFull(4:6,:),'+');
hold off
title('Root rotation position')

figure()
if model.nq == 12
subplot(211)
end
hold on
plot(t_estim,q_opt(7:9,:),'o');
plot(t_simu,data.xFull(7:9,:),'x');
plot(t_simu,data.gaussianNoiseXFull(7:9,:),'+');
hold off
title('Arm1 positions')

if model.nq == 12
subplot(212)
hold on
plot(t_estim,q_opt(10:12,:),'o');
plot(t_simu,data.xFull(10:12,:),'x');
plot(t_simu,data.gaussianNoiseXFull(10:12,:),'+');
hold off
title('Arm2 positions')
end

% VELOCITY PLOT
figure()
subplot(211)
hold on
plot(t_estim,v_opt(1:3,:),'o');
plot(t_simu,data.xFull(1+model.nq:3+model.nq,:),'x');
plot(t_simu,data.gaussianNoiseXFull(1+model.nq:3+model.nq,:),'+');
hold off
title('Root translation velocity')

subplot(212)
hold on
plot(t_estim,v_opt(4:6,:),'o');
plot(t_simu,data.xFull(4+model.nq:6+model.nq,:),'x');
plot(t_simu,data.gaussianNoiseXFull(4+model.nq:6+model.nq,:),'+');
hold off
title('Root rotation velocity')

figure()
if model.nq == 12
subplot(211)
end
hold on
plot(t_estim,v_opt(7:9,:),'o');
plot(t_simu,data.xFull(7+model.nq:9+model.nq,:),'x');
plot(t_simu,data.gaussianNoiseXFull(7+model.nq:9+model.nq,:),'+');
hold off
title('Arm1 velocity')

if model.nq == 12
subplot(212)
hold on
plot(t_estim,v_opt(10:12,:),'o');
plot(t_simu,data.xFull(10+model.nq:12+model.nq,:),'x');
plot(t_simu,data.gaussianNoiseXFull(10+model.nq:12+model.nq,:),'+');
hold off
title('Arm2 velocity')
end

% CONTROL PLOT
figure()
if model.nq == 12
subplot(211)
end
hold on
plot(t_estim,u_opt(1:3,:),'o');
plot(t_simu,data.uFull(1:3,:),'x');
hold off
title('Arms control')

if model.nq == 12
subplot(212)
hold on
plot(t_estim,u_opt(4:6,:),'o');
plot(t_simu,data.uFull(4:6,:),'x');
hold off
title('Arm2 control')
end
% showmotion(model, 0:data.Duration/(data.Nint-1):data.Duration, q_opt(:,:))

% L = @(q)base_referential_coor(model, q);
% S = @(u) 0.0003* (u'*u);
% 
% temp_markers = size(model.markers.coordinates);
% N_cardinal_coor = temp_markers(1);
% N_markers = temp_markers(2);
% markers = zeros(data.Nint, N_cardinal_coor * N_markers);
% 
% for old_value = 1:data.Nint
%     markers(old_value,:) = data.gaussianNoiseMarkers(old_value,:);
% end
% 
% J = full(sol.f);
% Jx = 0;
% Ju = 0;
% 
% T = data.Duration; % secondes
% N = data.Nint; % nb colloc nodes
% dN = T/N;
% M = 4;
% DT = dN/M;
% for i=1:data.Nint
%     for j=1:M
%         k1_q = S(u_opt(:,i));
%         k2_q = S(u_opt(:,i));
%         k3_q = S(u_opt(:,i));
%         k4_q = S(u_opt(:,i));
% 
%         Ju=Ju+DT/6*(k1_q +2*k2_q +2*k3_q +k4_q);
%         
%         k1_q = objective_func(model,data.gaussianNoiseMarkers(k,:),L(q_opt(:,i)));
%         k2_q = objective_func(model,data.gaussianNoiseMarkers(k,:),L(q_opt(:,i) + DT/M/2 * k1_q));
%         k3_q = objective_func(model,data.gaussianNoiseMarkers(k,:),L(q_opt(:,i) + DT/M/2 * k2_q));
%         k4_q = objective_func(model,data.gaussianNoiseMarkers(k,:),L(q_opt(:,i) + DT/M * k3_q));
%         
%         Jx=Jx+DT/6*(k1_q +2*k2_q +2*k3_q +k4_q);
%     end
% %     Ju = Ju + 0.5*u_opt(:,i)'*u_opt(:,i);
% end

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
% % end
% 
% function J = objective_func(model,markers,estimated_markers)
% J = 0;
% 
% temp_markers = size(model.markers.coordinates);
% N_cardinal_coor = temp_markers(1);
% N_markers = temp_markers(2);
% 
% n = 0;
% for m = 1:N_markers
%     distance_between_points = 0;
%     for l = 1:N_cardinal_coor
%         n = n + 1;
%         distance_between_points = distance_between_points + (markers(n) - estimated_markers{n}).^2;
%     end
%     J = J + 0.5 * distance_between_points;
% end
% end