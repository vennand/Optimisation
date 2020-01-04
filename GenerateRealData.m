function [model, data] = GenerateRealData(model,data)

% Assuming there is only one variable, and that it is q
real_data = load('Do_822_contact_2_MOD200.00_GenderF_DoCig_Q.mat');
real_data = struct2cell(real_data);
q = real_data{1};

Nint = data.Nint;
realNint = size(q);

for old_value = 1:Nint
    new_value = range_conversion(old_value, Nint, 1, realNint, 1);
    
    something = base_referential_coor(model,q(:,new_value)); % FINISH THIIIIIIS!!!
end




GaussianNoise_PosMarkers = zeros(simNint, N_cardinal_coor * N_markers);
for i = 1:simNint
    GaussianNoise_PosMarkers(i,:) = base_referential_coor(model,GaussianNoise_x(:,i));
end

% Storing data
data.xFull = Xi;
data.uFull = Ui;
data.markersFull = PosMarkers;
data.gaussianNoiseXFull = GaussianNoise_x;
data.gaussianNoiseMarkersFull = GaussianNoise_PosMarkers;

% Scaling data
Xi = zeros(model.nx,Nint);
Ui = zeros(model.nu,Nint);
PosMarkers = zeros(Nint, N_cardinal_coor * N_markers);
GaussianNoise_x = zeros(model.nx,Nint);
GaussianNoise_PosMarkers = zeros(Nint, N_cardinal_coor * N_markers);
for old_value = 1:Nint
    new_value = range_conversion(old_value, Nint, 1, simNint, 1);
    
    Xi(:,old_value) = data.xFull(:,floor(new_value));
    Ui(:,old_value) = data.uFull(:,floor(new_value));
    PosMarkers(old_value,:) = data.markersFull(floor(new_value),:);
    GaussianNoise_x(:,old_value) = data.gaussianNoiseXFull(:,floor(new_value));
    GaussianNoise_PosMarkers(old_value,:) = data.gaussianNoiseMarkersFull(floor(new_value),:);
end

data.x = Xi;
data.u = Ui;
data.markers = PosMarkers;
data.gaussianNoiseX = GaussianNoise_x;
data.gaussianNoiseMarkers = GaussianNoise_PosMarkers;
end

function new_value = range_conversion(old_value, old_max, old_min, new_max, new_min)
old_range = (old_max - old_min);
new_range = (new_max - new_min);
new_value = (((old_value - old_min) * new_range) / old_range) + new_min;
end