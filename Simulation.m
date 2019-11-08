% clear; close all; clc;

run('startup.m')
import casadi.*

model.NB = 4;
model.parent = [0 1 2 3];
model.jtype = {'R', 'Rx', 'Ry', 'Rz'};

model.Xtree = { eye(6), pluho([0 1 0 0;-1 0 0 1;0 0 1 0;0 0 0 1]), eye(6), eye(6) };

rod1 = mcI( 1, [0.5,0,0], diag([0.01,1,1]) );
rod2 = mcI( 1, [0.5,0,0], diag([0.01,1,1]) );
model.I = { rod1, zeros(6,6), zeros(6,6), rod2 };

model.markers.name = {'A','B','C','D','E','F'};
model.markers.parent = [6 6 6 9 9 9];
model.markers.coordinates = [1/4 0.05 0;2/4 0 -0.05;3/4 -0.05 0.05;...
                             1/4 0.05 0;2/4 0 -0.05;3/4 -0.05 0.05]';

model.appearance.base = { 'line', [1.1 0 0; 0 0 0; 0 1.1 0; 0 0 0; 0 0 1.1]};
model.appearance.body{1} = { %'cyl', [0 0 0; 1 0 0], 0.01, ...
                             'sphere', model.markers.coordinates(:,1),0.03, ...
                             'sphere', model.markers.coordinates(:,2),0.03, ...
                             'sphere', model.markers.coordinates(:,3),0.03};
model.appearance.body{2} = {};
model.appearance.body{3} = {};
model.appearance.body{4} = { %'cyl', [0 0 0; 1 0 0], 0.01, ...
                             'sphere', model.markers.coordinates(:,4),0.03, ...
                             'sphere', model.markers.coordinates(:,5),0.03, ...
                             'sphere', model.markers.coordinates(:,6),0.03};

model = floatbase(model);

% model.gravity = [0 0 0];

T = 1; % secondes
dt = 0.01; % secondes
nb_dt = T/dt;
N = model.NB;

q = zeros(N,1);
qd = zeros(N,1);
tau = zeros(N,1);

qd(1,1) = 1;
qd(3,1) = 9.81/2;
% qd(4,1) = 0;
qd(7,1) = 2;
qd(8,1) = 2;
qd(9,1) = 2;

f_ext = {};

dyn = zeros(nb_dt,N,3);

Xa = cell(N,1);
TransMatrix = cell(nb_dt,N);

N_markers = size(model.markers.coordinates);
PosMarkers = cell(nb_dt,N_markers(1), N_markers(2));

for i = 1:nb_dt
    if i == 1
        qdd = FDab( model, q, qd, tau, f_ext );
    else
        qd = dyn(i-1,:,2) + dyn(i-1,:,3)*dt;
        q = dyn(i-1,:,1) + qd*dt;
        qdd = FDab( model, q, qd, tau, f_ext );
    end
    
    dyn(i,:,1) = q;
    dyn(i,:,2) = qd;
    dyn(i,:,3) = qdd;
    
    PosMarkers(i,:,:) = base_referential_coor(model,q) ;
end

variance = 0.01;

gaussian_noise_dyn = dyn(:,:,1) + sqrt(variance).*(2.*rand(nb_dt,N)-1);
% plot(dt:dt:T,[gaussian_noise_dyn(:,3,1) dyn(:,3,1)])


GaussianNoise_PosMarkers = cell(nb_dt,N_markers(1), N_markers(2));
for i = 1:nb_dt
    GaussianNoise_PosMarkers(i,:,:) = base_referential_coor(model,gaussian_noise_dyn(i,:)) ;
end

% variance = 0.01;
% 
% GaussianNoise_PosMarkers = cell(nb_dt,N_markers(1), N_markers(2));
% for k = 1:N_markers(2)
%     noise = sqrt(variance).*(2.*rand(1,3)-1);
%     for i = 1:nb_dt
%         for j = 1:3
%             GaussianNoise_PosMarkers{i,j,k} = PosMarkers{i,j,k} + noise(j);
%         end
%     end
% end

% scatter3(dyn(:,1,1), dyn(:,2,1), dyn(:,3,1), 1)
% hold on
% scatter3(gaussian_noise_dyn(:,1,1), gaussian_noise_dyn(:,2,1), gaussian_noise_dyn(:,3,1), 1)

% showmotion(model, dt:dt:T, dyn(:,:,1).')
% showmotion(model)