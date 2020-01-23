fps = 0.1;%seconde


ise = evalin( 'base', 'exist(''q_opt'',''var'') == 1' );
if ise
    sol_markers = zeros(data.Nint+1,285);
    for i=1:data.Nint+1
        sol_markers(i,:) = base_referential_coor(model, q_opt(:,i));
    end
    markers_estim = sol_markers;
    disp('Estimation solution')
else
    kalman_markers = zeros(data.Nint+1,285);
    kalman_q = data.kalman_q;
    for i=1:data.Nint+1
        kalman_markers(i,:) = base_referential_coor(model, kalman_q(:,i));
    end
    markers_estim = kalman_markers;
    disp('Kalman solution')
end

markers_c3d = data.markers;

for i=1:data.Nint+1
    scatter3(markers_estim(i,1:3:end),markers_estim(i,2:3:end),markers_estim(i,3:3:end))
    hold on
    scatter3(markers_c3d(i,1:3:end),markers_c3d(i,2:3:end),markers_c3d(i,3:3:end))
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