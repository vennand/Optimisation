clear; close all; clc;

disp('Running simulation...')
run('Simulation.m')

run('startup.m')
import casadi.*

opti = casadi.Opti();

% T = MX.sym('T', 1); % seconds
T = 1; % seconds
N = 21*T; % nb colloc nodes

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

for i = 1:N
    opti.subject_to(abs(u(7:9,i))<ones(3,1)*10);% constraint on position
end

% Trapezoidal Collocation
% q = x(1:model.NB,:);
% qd = x(model.NB+1:end,:);
% tau = u;
f_ext = {};
fx = [];

temp_markers = size(model.markers.coordinates);
N_cardinal_coor = temp_markers(1);
N_markers = temp_markers(2);
new_range_posmarkers = cell(N, N_cardinal_coor, N_markers);
Estimated_PosMarkers = cell(N, N_cardinal_coor, N_markers);

for old_value = 1:N
    new_value = range_conversion(old_value, N, 1, nb_dt, 1);
    new_range_posmarkers(old_value,:,:) = GaussianNoise_PosMarkers(floor(new_value),:,:);
end

disp('Calculating dynamics...')
qdd = FDab( model, x(1:model.NB,1), x(model.NB+1:end,1), u(:,1), f_ext );
for n = 1:N-1
    fx = [x(model.NB+1:end,n);qdd];
    
    qdd = FDab( model,x(1:model.NB,n+1) , x(model.NB+1:end,n+1), u(:,n+1), f_ext );
    fx1 = [x(model.NB+1:end,n+1);qdd];
    
    xkend = x(:,n) + 1/2.*dN*(fx1+fx);

    opti.subject_to(xkend - x(:,n+1) == 0);
end

disp('Calculating objective function...')
for n = 1:N
    Estimated_PosMarkers(n,:,:) = base_referential_coor(model, x(1:model.NB,n));
    
    for m = 1:N_markers
        distance_between_points = 0;
        for l = 1:N_cardinal_coor
            distance_between_points = distance_between_points + (new_range_posmarkers{n,l,m} - Estimated_PosMarkers{n,l,m}).^2;
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

function new_value = range_conversion(old_value, old_max, old_min, new_max, new_min)
    old_range = (old_max - old_min);
    new_range = (new_max - new_min);
    new_value = (((old_value - old_min) * new_range) / old_range) + new_min;
end