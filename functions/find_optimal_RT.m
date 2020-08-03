function RotoTrans_opt = find_optimal_RT(model,data)
import casadi.*
tic
[~, N_markers] = size(model.markers.coordinates);

real_data = ezc3dRead(data.dataFile);

labels = data.labels;
labels_name = real_data.parameters.POINT.LABELS.DATA(labels);
[~, order] = ismember(model.markers.name, labels_name);

markers = real_data.data.points(:,labels,1)/1000; %meter
markers = [markers(:, order, :); ones(1,95)];

model_markers = [base_referential_coor(model,zeros(1,42)); ones(1,95)];

RotoTrans = @(transX, transY, transZ, angleX, angleY, angleZ) ...
            [cos(angleZ)*cos(angleY),   cos(angleZ)*sin(angleY)*sin(angleX)-sin(angleZ)*cos(angleX),    cos(angleZ)*sin(angleY)*cos(angleX)+sin(angleZ)*sin(angleX),    transX;...
             sin(angleZ)*cos(angleY),   sin(angleZ)*sin(angleY)*sin(angleX)+cos(angleZ)*cos(angleX),    sin(angleZ)*sin(angleY)*cos(angleX)-cos(angleZ)*sin(angleX),    transY;...
             -sin(angleY),              cos(angleY)*sin(angleX),                                        cos(angleY)*cos(angleX),                                        transZ;...
             0,                         0,                                                              0,                                                              1];

% syms angleX;
% syms angleY;
% syms angleZ;
% 
% RotX = [1, 0, 0;...
%         0, cos(angleX), -sin(angleX);...
%         0, sin(angleX), cos(angleX)];
% 
% 
% RotY = [cos(angleY), 0, sin(angleY);...
%         0, 1, 0;...
%         -sin(angleY), 0, cos(angleY)];
% 
%     
% RotZ = [cos(angleZ), -sin(angleZ), 0;...
%         sin(angleZ), cos(angleZ), 0;...
%         0, 0, 1];
% 
% Roto = RotZ  * RotY * RotX

% Start with an empty NLP
J=0;
w={};
lbw = [];
ubw = [];
g={};
lbg = [];
ubg = [];

transX = SX.sym('transX', 1);
transY = SX.sym('transY', 1);
transZ = SX.sym('transZ', 1);

angleX = SX.sym('angleX', 1);
angleY = SX.sym('angleY', 1);
angleZ = SX.sym('angleZ', 1);

w = {w{:}, transX};
lbw = [lbw; -10];
ubw = [ubw;  10];
w = {w{:}, transY};
lbw = [lbw; -10];
ubw = [ubw;  10];
w = {w{:}, transZ};
lbw = [lbw; -10];
ubw = [ubw;  10];

w = {w{:}, angleX};
lbw = [lbw; -pi];
ubw = [ubw;  pi];
w = {w{:}, angleY};
lbw = [lbw; -pi];
ubw = [ubw;  pi];
w = {w{:}, angleZ};
lbw = [lbw; -pi];
ubw = [ubw;  pi];

RotoTrans_sym = RotoTrans(transX, transY, transZ, angleX, angleY, angleZ);

for i = 1:N_markers
    if ~isnan(markers(1,i))
        J = J + sum((model_markers(:,i) - RotoTrans_sym*markers(:,i)).^2);
    end
end

w = vertcat(w{:});
prob = struct('f', J, 'x', w);

options = struct;
options.ipopt.max_iter = 3000;
options.ipopt.print_level = 5;
% options.ipopt.linear_solver = 'ma57';

% options.ipopt.tol = 1e-5; % default: 1e-08
% options.ipopt.acceptable_tol = 1e-4; % default: 1e-06
% options.ipopt.constr_viol_tol = 0.001; % default: 0.0001
% options.ipopt.acceptable_constr_viol_tol = 0.1; % default: 0.01

% options.ipopt.hessian_approximation = 'limited-memory';

% disp('Generating Solver')
solver = nlpsol('solver', 'ipopt', prob, options);

w0 = zeros(6,1);
w0(4) = 0.0480;
w0(5) = -0.0657;
w0(6) = 1.5720;

sol = solver('x0', w0, 'lbx', lbw, 'ubx', ubw);

w_opt = full(sol.x);
transX_opt = w_opt(1);
transY_opt = w_opt(2);
transZ_opt = w_opt(3);
angleX_opt = w_opt(4);
angleY_opt = w_opt(5);
angleZ_opt = w_opt(6);

RotoTrans_opt = RotoTrans(transX_opt, transY_opt, transZ_opt, angleX_opt, angleY_opt, angleZ_opt);
toc
% Plot
% markers_opt = RotoTrans_opt*markers;
% markers_old = refential_matrix()'*markers(1:3,:);
% scatter3(markers_opt(1,:),markers_opt(2,:),markers_opt(3,:),'x')
% hold on
% scatter3(model_markers(1,:),model_markers(2,:),model_markers(3,:),'o')
% scatter3(markers(1,:),markers(2,:),markers(3,:),'.')
% scatter3(markers_old(1,:),markers_old(2,:),markers_old(3,:),'+')
% axis equal
% %     set(gca,'visible','off')
% grid off
% xlabel('x')
% ylabel('y')
% zlabel('z')
% hold off
end