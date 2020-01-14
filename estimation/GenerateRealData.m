function [model, data] = GenerateRealData(model,data)
frames = data.frames;
labels = data.labels;

real_data = ezc3dRead(data.dataFile);
frequency = real_data.header.points.frameRate;
data.Duration = length(frames)/frequency;

markers = real_data.data.points(:,labels,frames)/1000; %meter
labels_name = real_data.parameters.POINT.LABELS.DATA(labels);

Nint = data.Nint;
[num_dimension, num_label, num_frames] = size(markers);

% Repositioning origin

markers_reformat = zeros(num_dimension, num_label, Nint+1);

for old_value = 1:Nint+1
    new_value = range_conversion(old_value, Nint+1, 1, num_frames, 1);
    markers_reformat(:,:,old_value) = markers(:,:,round(new_value));
end

[dummy, order] = ismember(model.markers.name, labels_name);
markers_reformat = markers_reformat(:,order,:);

markers_reformat = reshape(markers_reformat, num_dimension*num_label, Nint+1)';

% Storing data
data.markers = markers_reformat;
end