function [model, data] = GenerateRealData(model,data)
frames = data.frames;
labels = data.labels;

real_data = ezc3dRead(data.dataFile);
frequency = real_data.header.points.frameRate;
data.Duration = length(frames)/frequency;

markers = real_data.data.points(:,labels,frames)/1000; %meter

Nint = data.Nint;
[num_dimension, num_label, num_frames] = size(markers);

% Repositioning origin

% root_marker = model.markers.coordinates(1:3,1);
% translation_matrix = eye(4,4);
% translation_matrix(1:3,4) = root_marker;
% 
% real_data_origin = translation_matrix\[markers(:,1,1);1];
% real_data_origin = real_data_origin(1:3);

markers_reformat = zeros(num_dimension, num_label, Nint+1);

for old_value = 1:Nint+1
    new_value = range_conversion(old_value, Nint+1, 1, num_frames, 1);
    markers_reformat(:,:,old_value) = markers(:,:,round(new_value));
    
%     for label = 1:num_label
%         markers_reformat(:,label,old_value) = markers_reformat(:,label,old_value) - real_data_origin;
%     end
end

markers_reformat = reshape(markers_reformat, num_dimension*num_label, Nint+1)';

% Storing data
data.markers = markers_reformat;



% kalman_markers = zeros(Nint+1,285);
% kalman_q = data.kalman_q;
% for i=1:Nint+1
%     kalman_q(1:3,i) = [0;0;0];
%     kalman_markers(i,:) = base_referential_coor(model, kalman_q(:,i));
% end

% x = markers_reformat(1,1:3)';
% y = markers_reformat(1,4:6)';
% z = markers_reformat(1,7:9)';
% X = [x y z];
% 
% x_prime = kalman_markers(1,1:3)';
% y_prime = kalman_markers(1,4:6)';
% z_prime = kalman_markers(1,7:9)';
% X_prime = [x_prime y_prime z_prime];
% 
% R = X_prime / X;



% scatter3(kalman_markers(1,1:3:end),kalman_markers(1,2:3:end),kalman_markers(1,3:3:end))
% hold on
% scatter3(data.markers(1,1:3:end),data.markers(1,2:3:end),data.markers(1,3:3:end))
% scatter3(kalman_markers(5,1:3:end),kalman_markers(5,2:3:end),kalman_markers(5,3:3:end))
% scatter3(data.markers(5,1:3:end),data.markers(5,2:3:end),data.markers(5,3:3:end))
% scatter3(kalman_markers(10,1:3:end),kalman_markers(10,2:3:end),kalman_markers(10,3:3:end))
% scatter3(data.markers(10,1:3:end),data.markers(10,2:3:end),data.markers(10,3:3:end))
end