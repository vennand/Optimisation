function [] = AnimatePlot(model, data, cloud1, cloud2, fps)

is_fps = false;
if nargin == 5
    is_fps = true;
end

[N_cardinal_coor, N_markers] = size(model.markers.coordinates);

disp('Shapes are:')

switch cloud1
    case 'sol'
        sol_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
        for i=1:data.Nint+1
            sol_markers(:,:,i) = base_referential_coor(model, data.q_opt(:,i));
        end
        markers_cloud1 = sol_markers;
        disp(['Estimated cloud' blanks(4) 'o'])
    case 'kalman'
        kalman_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
        kalman_q = data.kalman_q;
        for i=1:data.Nint+1
            kalman_markers(:,:,i) = base_referential_coor(model, kalman_q(:,i));
        end
        markers_cloud1 = kalman_markers;
        disp(['Kalman cloud' blanks(4) 'o'])
    case 'kalmanUnopti'
        kalman_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
        kalman_q = data.kalman_qUnoptimised;
        for i=1:data.Nint+1
            kalman_markers(:,:,i) = base_referential_coor(model, kalman_q(:,i));
        end
        markers_cloud1 = kalman_markers;
        disp(['Kalman unoptimised cloud' blanks(4) 'o'])
    case 'mocap'
        markers_cloud1 = data.markers;
        disp(['Mocap cloud' blanks(4) 'o'])
    case 'sim'
        markers_cloud1 = data.markers_sim;
        disp(['Simulation cloud' blanks(4) 'o'])
    otherwise
        disp('No such dataset, n00b!')
        disp('Options are: sol, kalman, kalmanUnopti, mocap, sim')
end

switch cloud2
    case 'sol'
        sol_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
        for i=1:data.Nint+1
            sol_markers(:,:,i) = base_referential_coor(model, data.q_opt(:,i));
        end
        markers_cloud2 = sol_markers;
        disp(['Estimated cloud' blanks(4) 'x'])
    case 'kalman'
        kalman_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
        kalman_q = data.kalman_q;
        for i=1:data.Nint+1
            kalman_markers(:,:,i) = base_referential_coor(model, kalman_q(:,i));
        end
        markers_cloud2 = kalman_markers;
        disp(['Kalman cloud' blanks(4) 'x'])
    case 'kalmanUnopti'
        kalman_markers = zeros(N_cardinal_coor,N_markers,data.Nint+1);
        kalman_q = data.kalman_qUnoptimised;
        for i=1:data.Nint+1
            kalman_markers(:,:,i) = base_referential_coor(model, kalman_q(:,i));
        end
        markers_cloud2 = kalman_markers;
        disp(['Kalman unoptimised cloud' blanks(4) 'o'])
    case 'mocap'
        markers_cloud2 = data.markers;
        disp(['Mocap cloud' blanks(4) 'x'])
    case 'sim'
        markers_cloud2 = data.markers_sim;
        disp(['Simulation cloud' blanks(4) 'x'])
    otherwise
        disp('No such dataset, n00b!')
        disp('Options are: sol, kalman, kalmanUnopti, mocap, sim')
end

figure(20)
for i=1:data.Nint+1
    disp(['Node: ' num2str(i)])
    scatter3(markers_cloud1(1,:,i),markers_cloud1(2,:,i),markers_cloud1(3,:,i),'o')
    hold on
    scatter3(markers_cloud2(1,:,i),markers_cloud2(2,:,i),markers_cloud2(3,:,i),'x')
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

end