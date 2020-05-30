clear, clc, close all

load('Solutions/Do_822_F3100-3300_U1e-07_N50_optimisedKalman=1_optimisedKalmanGravity1_IPOPTMA57_.mat')
data_kalman_optim_no_grav = load(['Solutions/Do_822_F' num2str(data.frames(1)) '-' num2str(data.frames(end)) ...
                                  '_U' num2str(data.weightU) '_N' num2str(data.Nint) ...
                                  '_weightQV' num2str(1) '-' num2str(0.01) ...
                                  '_optimiseGravity=' num2str(false) ...
                                  '_gravityRotationBound=' num2str(data.gravityRotationBound) ...
                                  '_IPOPTMA57_Q.mat']);
data_kalman_optim_no_grav = data_kalman_optim_no_grav.('data');

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

J_estim = 0;
J_kalman_optimised = 0;
J_kalman_optimised_no_grav = 0;
J_kalman_unoptimised = 0;

distances_estim = NaN(N_cardinal_coor, N_markers,data.Nint+1);
distances_kalman_optimised = NaN(N_cardinal_coor, N_markers,data.Nint+1);
distances_kalman_optimised_no_grav = NaN(N_cardinal_coor, N_markers,data.Nint+1);
distances_kalman_unoptimised = NaN(N_cardinal_coor, N_markers,data.Nint+1);

estim_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
estim_q = data.q_opt;

kalman_markers_optimised = zeros(N_cardinal_coor,N_markers,data.Nint+1);
kalman_q_optimised = data.kalman_q;

kalman_markers_optimised_no_grav = zeros(N_cardinal_coor,N_markers,data.Nint+1);
kalman_q_optimised_no_grav = data_kalman_optim_no_grav.q_opt;

kalman_markers_unoptimised = zeros(N_cardinal_coor,N_markers,data.Nint+1);
kalman_q_unoptimised = data.kalman_qUnoptimised;

for i=1:data.Nint+1
    estim_markers(:,:,i) = base_referential_coor(model, estim_q(:,i));
    kalman_markers_optimised(:,:,i) = base_referential_coor(model, kalman_q_optimised(:,i));
    kalman_markers_optimised_no_grav(:,:,i) = base_referential_coor(model, kalman_q_optimised_no_grav(:,i));
    kalman_markers_unoptimised(:,:,i) = base_referential_coor(model, kalman_q_unoptimised(:,i));
end

markers = data.markers;
is_nan = double(isnan(markers));

for i=1:data.Nint+1
    for m = 1:N_markers
        distance_between_points_estim = 0;
        distance_between_points_kalman_optimised = 0;
        distance_between_points_kalman_optimised_no_grav = 0;
        distance_between_points_kalman_unoptimised = 0;
        for l = 1:N_cardinal_coor
            if ~is_nan(l,m,i)
                distance_between_points_estim = distance_between_points_estim + (markers(l,m,i) - estim_markers(l,m,i)).^2;
                distance_between_points_kalman_optimised = distance_between_points_kalman_optimised + (markers(l,m,i) - kalman_markers_optimised(l,m,i)).^2;
                distance_between_points_kalman_optimised_no_grav = distance_between_points_kalman_optimised_no_grav + (markers(l,m,i) - kalman_markers_optimised_no_grav(l,m,i)).^2;
                distance_between_points_kalman_unoptimised = distance_between_points_kalman_unoptimised + (markers(l,m,i) - kalman_markers_unoptimised(l,m,i)).^2;
            end
        end
        distances_estim(l,m,i) = 0.5 * sqrt(distance_between_points_estim);
        distances_kalman_optimised(l,m,i) = 0.5 * sqrt(distance_between_points_kalman_optimised);
        distances_kalman_optimised_no_grav(l,m,i) = 0.5 * sqrt(distance_between_points_kalman_optimised_no_grav);
        distances_kalman_unoptimised(l,m,i) = 0.5 * sqrt(distance_between_points_kalman_unoptimised);
        
        J_estim = J_estim + 0.5 * distance_between_points_estim;
        J_kalman_optimised = J_kalman_optimised + 0.5 * distance_between_points_kalman_optimised;
        J_kalman_optimised_no_grav = J_kalman_optimised_no_grav + 0.5 * distance_between_points_kalman_optimised_no_grav;
        J_kalman_unoptimised = J_kalman_unoptimised + 0.5 * distance_between_points_kalman_unoptimised;
    end
end

disp('Estimation with optimised Kalman')
disp(nansum(distances_estim,'all')/N_markers/data.Nint)
disp('Kalman optimised')
disp(nansum(distances_kalman_optimised,'all')/N_markers/data.Nint)
disp('Kalman optimised without gravity')
disp(nansum(distances_kalman_optimised_no_grav,'all')/N_markers/data.Nint)
disp('Kalman unoptimised')
disp(nansum(distances_kalman_unoptimised,'all')/N_markers/data.Nint)