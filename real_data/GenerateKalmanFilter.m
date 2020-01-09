function [model, data] = GenerateKalmanFilter(model,data)
frames = data.frames;

% Assuming there is only one variable, and that it is q
kalmanData = load(data.kalmanDataFile);
kalmanData = struct2cell(kalmanData);
q = kalmanData{1};
q = q(:,frames);

Nint = data.Nint;
realNint = data.realNint;
% [N_cardinal_coor, N_markers] = size(model.markers.coordinates);
new_q = zeros(dof, Nint);
% PosMarkers = zeros(Nint, N_cardinal_coor * N_markers);

for old_value = 1:Nint
    new_value = range_conversion(old_value, Nint, 1, realNint, 1);
    new_q(:,old_value) = q(:,floor(new_value));
end

% for old_value = 1:Nint
%     new_value = range_conversion(old_value, Nint, 1, realNint, 1);
%     PosMarkers(old_value,:) = base_referential_coor(model,q(:,floor(new_value)));
% end

% Storing data
data.kalman_q = new_q;
data.kalman_qFull = q;
% data.markers = PosMarkers;

end