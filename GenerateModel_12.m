% Only works with nDoF = 12
function [model, data] = GenerateModel_12(nDoF,data)
import casadi.*

nDoF = '12';

model.name = nDoF;

model.NB = str2double(nDoF)-5; % model.floatbase will add 5 degrees of freedom
model.parent = [0 1 2 3 1 5 6];
model.jtype = {'R', 'Rx', 'Ry', 'Rz', 'Rx', 'Ry', 'Rz'};

model.Xtree = { eye(6), pluho([0 1 0 0;-1 0 0 1;0 0 1 0;0 0 0 1]), eye(6), eye(6), pluho([0 1 0 0;-1 0 0 0;0 0 1 0;0 0 0 1]), eye(6), eye(6) };

rod1 = mcI( 1, [0.5,0,0], diag([0.01,1,1]) );
rod2 = mcI( 1, [0.5,0,0], diag([0.01,1,1]) );
rod3 = mcI( 1, [0.5,0,0], diag([0.01,1,1]) );
model.I = { rod1, zeros(6,6), zeros(6,6), rod2, zeros(6,6), zeros(6,6), rod3  };

model.markers.name = {'A','B','C','D','E','F','G','H','I'};
model.markers.parent = [6 6 6 9 9 9 12 12 12 ]; % Yeah??
model.markers.coordinates = [1/4 0.05 0;2/4 0 -0.05;3/4 -0.05 0.05;...
                             1/4 0.05 0;2/4 0 -0.05;3/4 -0.05 0.05;...
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
model.appearance.body{5} = {};
model.appearance.body{6} = {};
model.appearance.body{7} = { %'cyl', [0 0 0; 1 0 0], 0.01, ...
                             'sphere', model.markers.coordinates(:,7),0.03, ...
                             'sphere', model.markers.coordinates(:,8),0.03, ...
                             'sphere', model.markers.coordinates(:,9),0.03};

model = floatbase(model);

model.nq = model.NB;
model.nx = model.nq+model.nq;
model.nu = model.nq-6;

data.x0 = zeros(model.nx,1);
data.u0 = zeros(model.nu,1);

model.idx_q = 1:model.nq;
model.idx_v = model.nq+1:2*model.nq;

xmin_base = [-inf,-inf,-inf,-inf,-inf,-inf];
xmax_base = [ inf, inf, inf, inf, inf, inf];

model.xmin   =  [xmin_base'; -inf*ones(model.nq-6,1); ...
                 xmin_base'; -100*ones(model.nq-6,1)];
model.xmax   =  [xmax_base';  inf*ones(model.nq-6,1); ...
                 xmax_base';  100*ones(model.nq-6,1)];

model.umin = -50*ones(model.nu,1);
model.umax =  50*ones(model.nu,1);


% model.gravity = [0 0 0];

end