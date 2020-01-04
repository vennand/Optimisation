function [model, data] = GenerateSimulation_RK4(model,data)
import casadi.*

simNint = data.simNint;
Nint = data.Nint;

M = 4;
T = data.Duration; % secondes
DT = T/simNint/M; % secondes

N = model.NB;

x = SX.sym('x', model.nx,1);
u = SX.sym('u', model.nu,1);

tau_base = SX.zeros(6,1);
forDyn = @(x,u)[  x(model.idx_v)
    FDab_Casadi( model, x(model.idx_q), x(model.idx_v), vertcat(tau_base,u) )];

f = Function('f', {x, u}, {forDyn(x,u)});

x = zeros(model.nx,1);
u = zeros(model.nu,1);

% Initial velocities
x(N+1) = 1;
x(N+3) = 9.81/2;
% x(N+4) = 0;
% x(N+7) = 2;
% x(N+8) = 2;
% x(N+9) = 2;

% u(1) = 0.2;

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);
PosMarkers = zeros(simNint, N_cardinal_coor * N_markers);

Xi = zeros(model.nx,simNint);
Ui = zeros(model.nu,simNint);
Xi(:,1) = x;
Ui(:,1) = u;
PosMarkers(1,:) = base_referential_coor(model,x(model.idx_q));
for i = 2:simNint
    for j=1:M
        k1 = f(x, u);
        k2 = f(x + DT/2 * k1, u);
        k3 = f(x + DT/2 * k2, u);
        k4 = f(x + DT * k3, u);
        
        x=x+DT/6*(k1 +2*k2 +2*k3 +k4);
    end
    x = full(x);
    Xi(:,i) = x;
    
%     if i == floor(simNint/8-1)
%         u(4) = -0.2;
%     elseif i == floor(simNint*3/8-1)
%         u(1) = 0.2;
%     elseif i == floor(simNint*4/8-1)
%         u(4) = 0.2;
%     elseif i == floor(simNint*5/8-1)
%         u(1) = -0.2;
%     elseif i == floor(simNint*6/8-1)
%         u(4) = -0.2;
%     elseif i == floor(simNint*7/8-1)
%         u(1) = 0.2;
%     end
    
    Ui(:,i) = u;
    
    PosMarkers(i,:) = base_referential_coor(model,x(model.idx_q));
end

variance = data.simVariance;
rng(0,'twister'); % So that the gaussian is always the same
GaussianNoise_x = Xi + sqrt(variance).*(2.*rand(model.nx,simNint)-1);

GaussianNoise_PosMarkers = zeros(simNint, N_cardinal_coor * N_markers);
for i = 1:simNint
    GaussianNoise_PosMarkers(i,:) = base_referential_coor(model,GaussianNoise_x(:,i));
end

% Storing data
data.xFull = Xi;
data.uFull = Ui;
data.markersFull = PosMarkers;
data.gaussianNoiseXFull = GaussianNoise_x;
data.gaussianNoiseMarkersFull = GaussianNoise_PosMarkers;

% Scaling data
Xi = zeros(model.nx,Nint);
Ui = zeros(model.nu,Nint);
PosMarkers = zeros(Nint, N_cardinal_coor * N_markers);
GaussianNoise_x = zeros(model.nx,Nint);
GaussianNoise_PosMarkers = zeros(Nint, N_cardinal_coor * N_markers);
for old_value = 1:Nint
    new_value = range_conversion(old_value, Nint, 1, simNint, 1);
    
    Xi(:,old_value) = data.xFull(:,floor(new_value));
    Ui(:,old_value) = data.uFull(:,floor(new_value));
    PosMarkers(old_value,:) = data.markersFull(floor(new_value),:);
    GaussianNoise_x(:,old_value) = data.gaussianNoiseXFull(:,floor(new_value));
    GaussianNoise_PosMarkers(old_value,:) = data.gaussianNoiseMarkersFull(floor(new_value),:);
end

data.x = Xi;
data.u = Ui;
data.markers = PosMarkers;
data.gaussianNoiseX = GaussianNoise_x;
data.gaussianNoiseMarkers = GaussianNoise_PosMarkers;
end

function new_value = range_conversion(old_value, old_max, old_min, new_max, new_min)
old_range = (old_max - old_min);
new_range = (new_max - new_min);
new_value = (((old_value - old_min) * new_range) / old_range) + new_min;
end