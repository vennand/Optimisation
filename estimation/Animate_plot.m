fps = 1;%seconde

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

ise = evalin( 'base', 'exist(''q_opt'',''var'') == 1' );
if ise
    sol_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
    for i=1:data.Nint+1
        sol_markers(:,:,i) = base_referential_coor(model, q_opt(:,i));
    end
    markers_estim = sol_markers;
    disp('Estimation solution')
else
    kalman_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
    kalman_q = data.kalman_q;
    for i=1:data.Nint+1
        kalman_markers(:,:,i) = base_referential_coor(model, kalman_q(:,i));
    end
    markers_estim = kalman_markers;
    disp('Kalman solution')
end

markers_c3d = data.markers;

for i=1:data.Nint+1
    scatter3(markers_estim(1,:,i),markers_estim(2,:,i),markers_estim(3,:,i))
    hold on
    scatter3(markers_c3d(1,:,i),markers_c3d(2,:,i),markers_c3d(3,:,i))
    axis equal
%     set(gca,'visible','off')
    grid off
    xlabel('x')
    ylabel('y')
    zlabel('z')
    drawnow
    pause(fps)
    hold off
end