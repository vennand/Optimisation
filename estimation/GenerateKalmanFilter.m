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

if data.optimisedKalman
    optimisedKalman_filename = ['Solutions/Do_822_F' num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
                                '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
                                '_weightQV' num2str(1) '-' num2str(0.01) ...
                                '_optimiseGravity=' num2str(data.optimisedKalmanGravity) ...
                                '_gravityRotationBound=' num2str(data.gravityRotationBound) ...
                                '_IPOPTMA57_Q.mat'];
    if isfile(optimisedKalman_filename)
        data_kalman = load(optimisedKalman_filename, 'data');
        data_kalman = data_kalman.('data');
    else
        error('No EKF found O.o')
    end
    
    % Storing data
    data.kalman_q = data_kalman.q_opt;
    data.kalman_v = data_kalman.v_opt;
    data.kalman_tau = data_kalman.u_opt;
    
    data.kalman_qUnoptimised = new_q;
    data.kalman_vUnoptimised = new_v;
    data.kalman_aUnoptimised = new_a;
    data.kalman_tauUnoptimised = new_tau(7:end,1:end-1);
    
    if data.optimisedKalmanGravity
        data.gravity = data_kalman.G_opt;
    end
else
    data.kalman_q = new_q;
    data.kalman_v = new_v;
    data.kalman_a = new_a;
    data.kalman_tau = new_tau(7:end,1:end-1);
end

data.kalman_qFull = q;
data.kalman_vFull = v;
data.kalman_aFull = a;
data.kalman_tauFull = tau(:,1:end-1);

end