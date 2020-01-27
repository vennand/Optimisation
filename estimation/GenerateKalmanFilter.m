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

% Nint = data.Nint;
% realNint = data.realNint;
% 
% new_q = zeros(dof, Nint+1);
% new_v = zeros(dof, Nint+1);
% new_a = zeros(dof, Nint+1);
% new_tau = zeros(dof-6, Nint+1);
% 
% for old_value = 1:Nint+1
%     new_value = range_conversion(old_value, Nint+1, 1, realNint, 1);
%     
%     new_q(:,old_value) = q(:,round(new_value));
%     new_v(:,old_value) = v(:,round(new_value));
%     new_a(:,old_value) = a(:,round(new_value));
%     new_tau(:,old_value) = tau(7:end,round(new_value));
% end

new_q = q(:, 1:data.step:end);
new_v = v(:, 1:data.step:end);
new_a = a(:, 1:data.step:end);
new_tau = tau(:, 1:data.step:end);

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