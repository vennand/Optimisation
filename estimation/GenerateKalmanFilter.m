function [model, data] = GenerateKalmanFilter(model,data)
frames = data.frames;

% Assuming there is only one variable per file, and that they are q, v and a
kalmanData_q = load(data.kalmanDataFile_q);
kalmanData_v = load(data.kalmanDataFile_v);
kalmanData_a = load(data.kalmanDataFile_a);

kalmanData_q = struct2cell(kalmanData_q);
kalmanData_v = struct2cell(kalmanData_v);
kalmanData_a = struct2cell(kalmanData_a);

q = kalmanData_q{1};
v = kalmanData_v{1};
a = kalmanData_a{1};

q = q(:,frames);
v = v(:,frames);
a = a(:,frames);

[dof, frame_length] = size(q);

tau = zeros(dof,frame_length);
for i = 1:frame_length
    tau(:,i) = ID(model, q(:,i), v(:,i), a(:,i));
end

new_q = q(:, 1:data.step:end);
new_v = v(:, 1:data.step:end);
new_a = a(:, 1:data.step:end);
new_tau = tau(:, 1:data.step:end);

% Swap rotation of arms from xyz to yxz, to correct difference in models
% TEMPORARY!!!!
% new_q(15:16, :) = new_q(16:-1:15, :);
% new_q(24:25, :) = new_q(25:-1:24, :);
% new_v(15:16, :) = new_v(16:-1:15, :);
% new_v(24:25, :) = new_v(25:-1:24, :);
% new_a(15:16, :) = new_a(16:-1:15, :);
% new_a(24:25, :) = new_a(25:-1:24, :);
% new_tau(15:16, :) = new_tau(16:-1:15, :);
% new_tau(24:25, :) = new_tau(25:-1:24, :);

% Storing data
data.kalman_q = new_q;
data.kalman_v = new_v;
data.kalman_a = new_a;
data.kalman_tau = new_tau(7:end,1:end-1);

data.kalman_qFull = q;
data.kalman_vFull = v;
data.kalman_aFull = a;
data.kalman_tauFull = tau(:,1:end-1);

end