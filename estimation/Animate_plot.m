run('../startup.m')

% fps = 0.1;%seconde

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

is_q_opt = evalin( 'base', 'exist(''q_opt'',''var'') == 1' );
is_fps = evalin( 'base', 'exist(''fps'',''var'') == 1' );

if is_q_opt
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
    disp(['Node: ' num2str(i)])
    scatter3(markers_estim(1,:,i),markers_estim(2,:,i),markers_estim(3,:,i),'o')
    hold on
    scatter3(markers_c3d(1,:,i),markers_c3d(2,:,i),markers_c3d(3,:,i),'x')
    axis equal
%     set(gca,'visible','off')
    grid off
    xlabel('x')
    ylabel('y')
    zlabel('z')
    title(['Node: ' num2str(i)])
    drawnow
    if is_fps
        pause(fps)
    else
        try
            while true
                k = waitforbuttonpress;
                if k == 1 % key stroke = 1, click = 0
                    value = double(get(gcf,'CurrentCharacter'));
                    switch value
                        case 13 % return key
                            break
                        case 32 % space key
                            break
                    end
                end
            end
        catch
            break
        end
    end
    hold off
end