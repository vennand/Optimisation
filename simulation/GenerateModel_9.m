% Only works with nDoF = 9
function [model, data] = GenerateModel(data)
import casadi.*

nDoF = '9';

model.name = nDoF;

model.NB = str2double(nDoF)-5;
model.parent = [0 1 2 3];
model.jtype = {'R', 'Rx', 'Ry', 'Rz'};

model.Xtree = { eye(6), inv(pluho([0 -1 0 1;1 0 0 0;0 0 1 0;0 0 0 1])), eye(6), eye(6) };

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

model.nq = model.NB;
model.nx = model.nq+model.nq;
model.nu = model.nq-6;

data.x0 = zeros(model.nx,1);
data.u0 = zeros(model.nu,1);

model.idx_q = 1:model.nq;
model.idx_v = model.nq+1:2*model.nq;

model.xmin   =  [-inf,-inf,-inf,-inf,-inf,-inf, -2*pi, -2*pi, -2*pi, ...
                 -inf,-inf,-inf,-inf,-inf,-inf, -100,  -100,  -100]';
model.xmax   =  [ inf, inf, inf, inf, inf, inf,  2*pi,  2*pi,  2*pi, ...
                  inf, inf, inf, inf, inf, inf,  100,   100,   100]';

model.umin = -50*ones(model.nu,1);
model.umax =  50*ones(model.nu,1);


% model.gravity = [0 0 0];

end
