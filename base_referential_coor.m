% From the coordinates of a point in a segment referential, finds the coordinates in the base referential

function PosMarkers = base_referential_coor(model, q)

    N = model.NB;

    Xa = cell(N,1);
    TransMatrix = cell(N,1);

    N_markers = size(model.markers.coordinates);
    PosMarkers = cell(N_markers(1), N_markers(2));
    
    % Get the rotational matrix to change from a segment referential to the
    % base referential
    for j = 1:N
        XJ = jcalc( model.jtype{j}, q(j) );
        Xa{j} = XJ * model.Xtree{j};
        if model.parent(j) ~= 0
          Xa{j} = Xa{j} * Xa{model.parent(j)};
        end
        Transform = pluho(Xa{j});
        TransMatrix{j} = inv(Transform);		% displacement is inverse of coord xform
    end

    for k = 1:N_markers(2)
        if model.markers.parent(k) == 6
            trackpos = [model.markers.coordinates(:,k)' 1] * TransMatrix{6}'; % Check for floating base?
        elseif model.markers.parent(k) == 9
            trackpos = [model.markers.coordinates(:,k)' 1] * TransMatrix{9}'; % Modify/improve this!
        end
        PosMarkers(:,k) = num2cell(trackpos(1:3));
    end
end