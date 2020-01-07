function [model, data] = GenerateRealData(model,data)

% Assuming there is only one variable, and that it is q
real_data = load(data.dataFile);
real_data = struct2cell(real_data);
q = real_data{1};

Nint = data.Nint;
realNint = size(q);
PosMarkers = zeros(Nint, N_cardinal_coor * N_markers);

for old_value = 1:Nint
    new_value = range_conversion(old_value, Nint, 1, realNint, 1);
    PosMarkers(old_value,:) = base_referential_coor(model,q(:,new_value)); % FINISH THIIIIIIS!!!
end

% Storing data
data.q = q;
data.markers = PosMarkers;

end

function new_value = range_conversion(old_value, old_max, old_min, new_max, new_min)
old_range = (old_max - old_min);
new_range = (new_max - new_min);
new_value = (((old_value - old_min) * new_range) / old_range) + new_min;
end