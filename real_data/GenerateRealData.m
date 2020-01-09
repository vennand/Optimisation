function [model, data] = GenerateRealData(model,data)
frames = data.frames;
labels = data.labels;

real_data = ezc3dRead(data.dataFile);
frequency = real_data.parameters.POINT.RATE.DATA;
data.Duration = length(frames)/frequency;

Nint = data.Nint;
[N_cardinal_coor, N_markers] = size(model.markers.coordinates);
PosMarkers = zeros(Nint, N_cardinal_coor * N_markers);

markers = real_data.data.points(:,labels,frames);
[num_dimension, num_label, num_frames] = size(markers);

% C'est vraiment dégueux ce code...
markers_reformat = zeros(num_frames, num_dimension * num_label);
marker = zeros(1,num_dimension);
for frame = 1:num_frames
    for label = 1:num_label
        for dimension = 1:num_dimension
            marker(dimension) = markers(dimension,label,frame);
        end
        section = (label-1)*num_dimension+1:(label-1)*num_dimension+num_dimension;
        markers_reformat(frame,section) = marker;
    end
end

for old_value = 1:Nint
    new_value = range_conversion(old_value, Nint, 1, num_frames, 1);
    PosMarkers(old_value,:) = markers_reformat(floor(new_value),:);
end

% Storing data
data.markers = PosMarkers;
end