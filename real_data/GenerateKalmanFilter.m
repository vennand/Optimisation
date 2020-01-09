function [model, data] = GenerateKalmanFilter(model,data)
frames = data.frames;

% Assuming there is only one variable, and that it is q
kalmanData = load(data.kalmanDataFile);
kalmanData = struct2cell(kalmanData);
q = kalmanData{1};
q = q(:,frames);

[dof, frame_length] = size(q);

Nint = data.Nint;
realNint = data.realNint;
new_q = zeros(dof, Nint);

for old_value = 1:Nint
    new_value = range_conversion(old_value, Nint, 1, realNint, 1);
    new_q(:,old_value) = q(:,floor(new_value));
end

% Storing data
data.kalman_q = new_q;
data.kalman_qFull = q;

end