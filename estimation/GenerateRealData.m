function [model, data] = GenerateRealData(model,data)
frames = data.frames;
labels = data.labels;

real_data = ezc3dRead(data.dataFile);
frequency = real_data.parameters.POINT.RATE.DATA;
data.Duration = length(frames)/frequency;

Nint = data.Nint;
[N_cardinal_coor, N_markers] = size(model.markers.coordinates);
PosMarkers = zeros(Nint+1, N_cardinal_coor * N_markers);

markers = real_data.data.points(:,labels,frames);
[num_dimension, num_label, num_frames] = size(markers);

markers_reformat = reshape(markers, num_dimension*num_label, num_frames)';

for old_value = 1:Nint+1
    new_value = range_conversion(old_value, Nint+1, 1, num_frames, 1);
    PosMarkers(old_value,:) = markers_reformat(round(new_value),:);
end

% Storing data
data.markers = PosMarkers;

% scatter3(markers_reformat(:,1),markers_reformat(:,2),markers_reformat(:,3),'.')
% hold on
% scatter3(PosMarkers(:,1),PosMarkers(:,2),PosMarkers(:,3),'o')
end