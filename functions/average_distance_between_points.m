% Calculate distance average distance between end of chain markers
clear, clc, close all
run('../startup.m')

opti_Q = load('/home/andre/Optimisation/gravity/Solutions/Do_822_F3100-3300_U1e-07_N50_weightQV1-0.01_gravityRotationBound=0.19635_IPOPTMA57_Q.mat');
opti_Q_EndChainMarkers = load('/home/andre/Optimisation/gravity/Solutions/Do_822_F3100-3300_U1e-07_N50_weightQV1-0.01_gravityRotationBound=0.19635_IPOPTMA57_Q_EndChainMarkers.mat');

opti_Q_data = opti_Q.data;
opti_Q_model = opti_Q.model;

opti_Q_EndChainMarkers_data = opti_Q_EndChainMarkers.data;
opti_Q_EndChainMarkers_model = opti_Q_EndChainMarkers.model;

N_cardinal_coor = opti_Q_data.nCardinalCoor;
marker_num = [37:39 59:61 76:78 93:95];
N_markers = length(marker_num);

kalman_markers = zeros(N_cardinal_coor,95,opti_Q_data.Nint+1);
for i=1:opti_Q_data.Nint+1
    kalman_markers(:,:,i) = base_referential_coor(opti_Q_model, opti_Q_data.kalman_q(:,i));
end

sol_markers_Q = zeros(N_cardinal_coor,95,opti_Q_data.Nint+1);
for i=1:opti_Q_data.Nint+1
    sol_markers_Q(:,:,i) = base_referential_coor(opti_Q_model, opti_Q_data.q_opt(:,i));
end

sol_markers_Q_EndChainMarkers = zeros(N_cardinal_coor,95,opti_Q_EndChainMarkers_data.Nint+1);
for i=1:opti_Q_EndChainMarkers_data.Nint+1
    sol_markers_Q_EndChainMarkers(:,:,i) = base_referential_coor(opti_Q_EndChainMarkers_model, opti_Q_EndChainMarkers_data.q_opt(:,i));
end

distance_between_kalman_markers = zeros(N_markers, opti_Q_data.Nint+1);
distance_between_sol_markers_Q = zeros(N_markers, opti_Q_data.Nint+1);
distance_between_sol_markers_Q_EndChainMarkers = zeros(N_markers, opti_Q_EndChainMarkers_data.Nint+1);

for i = 1:opti_Q_data.Nint+1
    for m = 1:N_markers
        point1 = kalman_markers(:,marker_num(m),i);
        point2 = opti_Q_data.markers(:,marker_num(m),i);
        point3 = sol_markers_Q(:,marker_num(m),i);
        point4 = sol_markers_Q_EndChainMarkers(:,marker_num(m),i);

        if isnan(point2(1)) || isnan(point3(1)) || isnan(point4(1))
            break
        end

        distance_between_kalman_markers(m,i) = sqrt(sum((point1 - point2).^2))*1000;
        distance_between_sol_markers_Q(m,i) = sqrt(sum((point3 - point2).^2))*1000;
        distance_between_sol_markers_Q_EndChainMarkers(m,i) = sqrt(sum((point4 - point2).^2))*1000;
    end
end

figure(1)
subplot(311)
plot(distance_between_kalman_markers')
title("Kalman vs. Mocap")
legend("Right hand X", "Right hand Y", "Right hand Z", "Left hand X", "Left hand Y", "Left hand Z", ...
       "Right foot X", "Right foot Y", "Right foot Z", "Left foot X", "Left foot Y", "Left foot Z")
   
subplot(312)
plot(distance_between_sol_markers_Q')
title("Optimisation Q vs. Mocap")

subplot(313)
plot(distance_between_sol_markers_Q_EndChainMarkers')
title("Optimisation Q et EndChainMarkers vs. Mocap")