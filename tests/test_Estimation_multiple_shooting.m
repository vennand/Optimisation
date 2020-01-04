clear; close all; clc;

disp('Running simulation...')
run('Simulation.m')

run('startup.m')
import casadi.*

opti = casadi.Opti();

% T = MX.sym('T', 1); % seconds
T = 1; % seconds
N = 10*T; % nb colloc nodes

dN = T/N;

x = opti.variable(2*model.NB,N);% state (q1,q2,qd1,qd2)
u = opti.variable(model.NB,N);% control (u)

J = 0; % initialization

opti.set_initial(x,zeros(2*model.NB,N))
% opti.set_initial(u,zeros(model.NB,N))
%opti.subject_to(x(:,1)==zeros(2*model.NB,1));% initial constraint
% opti.subject_to(x(10,1)==1);% initial constraint
% opti.subject_to(x(12,1)==9.81/2);% initial constraint
%opti.subject_to(u(:,1)==zeros(model.NB,1));% initial constraint
opti.subject_to(u(1:6,:)==zeros(6,N));% initial constraint

% final_position = [dyn(1000,:,1).';dyn(1000,:,2).'];
% opti.subject_to(x(:,end)==final_position);% terminal constraint on position

% opti.subject_to(abs(x(1,:))<5);% constraint on position

% Trapezoidal Collocation
% q = x(1:model.NB,:);
% qd = x(model.NB+1:end,:);
% tau = u;
f_ext = {};
fx = [];

temp_markers = size(model.markers.coordinates);
N_cardinal_coor = temp_markers(1);
N_markers = temp_markers(2);
Estimated_PosMarkers = cell(N, N_cardinal_coor, N_markers);

disp('Calculating dynamics...')
qdd = FDab( model, x(1:model.NB,1), x(model.NB+1:end,1), u(:,1), f_ext );
for n = 1:N-1
    x_e = rg4(model,x(:,n),u(:,n),N,T);

    opti.subject_to(x_e - x(:,n+1) == 0);
end

disp('Calculating objective function...')
for n = 1:N
    Estimated_PosMarkers(n,:,:) = base_referential_coor(model, x(1:model.NB,n));
    
    for m = 1:N_markers
        distance_between_points = 0;
        for l = 1:N_cardinal_coor
            distance_between_points = distance_between_points + (GaussianNoise_PosMarkers{ceil(n/N*nb_dt),l,m} - Estimated_PosMarkers{n,l,m}).^2;
        end
%         distance_between_points = sqrt(distance_between_points);
        J = J + 1/2 * distance_between_points;
    end
end

opti.minimize(J);% minimize objective function

opti.solver('ipopt');
sol = opti.solve();
solx = sol.value(x);
solu = sol.value(u);

t_opti = linspace(dN,T,N);
t_simu = linspace(dt,T,nb_dt);

% % POSITION PLOT
% figure()
% subplot(311)
% hold on
% plot(t_opti,solx(1:3,:),'o');
% plot(t_simu,dyn(:,1:3,1),'x');
% hold off
% title('Root translation position')
% 
% subplot(312)
% hold on
% plot(t_opti,solx(4:6,:),'o');
% plot(t_simu,dyn(:,4:6,1),'x');
% hold off
% title('Root rotation position')
% 
% subplot(313)
% hold on
% plot(t_opti,solx(7:9,:),'o');
% plot(t_simu,dyn(:,7:9,1),'x');
% hold off
% title('Arm positions')

% POSITION PLOT
figure()
subplot(311)
hold on
plot(t_opti,solx(1:3,:),'o');
plot(t_simu,gaussian_noise_dyn(:,1:3),'x');
hold off
title('Root translation position')

subplot(312)
hold on
plot(t_opti,solx(4:6,:),'o');
plot(t_simu,gaussian_noise_dyn(:,4:6),'x');
hold off
title('Root rotation position')

subplot(313)
hold on
plot(t_opti,solx(7:9,:),'o');
plot(t_simu,gaussian_noise_dyn(:,7:9),'x');
hold off
title('Arm positions')

% VELOCITY PLOT
figure()
subplot(311)
hold on
plot(t_opti,solx(10:12,:),'o');
hold off
title('Root translation velocity')

subplot(312)
hold on
plot(t_opti,solx(13:15,:),'o');
hold off
title('Root rotation velocity')

subplot(313)
hold on
plot(t_opti,solx(16:18,:),'o');
hold off
title('Arm velocity')

% CONTROL PLOT
figure()
subplot(311)
plot(t_opti,solu(1:3,:),'o');
title('Root translation control')

subplot(312)
plot(t_opti,solu(4:6,:),'o');
title('Root rotation control')

subplot(313)
plot(t_opti,solu(7:9,:),'o');
title('Arms control')

% showmotion(model, dN:dN:T, solx(1:9,:))

function x = rg4(model,x,u,N,T)
    import casadi.*
    
    forDyn = @(x,u)[  x(model.NB+1:end)
        FDab( model, x(1:model.NB), x(model.NB+1:end), u )];
    
    M = 4; % RK4 steps per interval
    DT = T/(N-1)/M;
    for j=1:M
       k1 = forDyn(x,u);
       k2 = forDyn(x + DT/M/2 * k1,u);
       k3 = forDyn(x + DT/M/2 * k2,u);
       k4 = forDyn(x + DT/M * k3,u);
       x=x+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
end
