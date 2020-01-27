% From the coordinates of a point in a segment referential, finds the coordinates in the base referential
function PosMarkers = base_referential_coor(model, q)
    import casadi.*
    
    N = model.NB;
    
    isSX = false;
    for i = 1:N
        if class(model.Xtree{i}) == "casadi.SX"
            isSX = true;
        end
    end

    isMX = false;

    [N_cardinal_coor, N_markers] = size(model.markers.coordinates);
    
    if class(q) == "casadi.SX" || isSX
        PosMarkers = SX.sym('marker', N_cardinal_coor, N_markers);
    elseif class(q) == "casadi.MX"
        PosMarkers = {};
        isMX = true;
    elseif class(q) == "casadi.DM"
        PosMarkers = DM(N_cardinal_coor, N_markers);
    else
        PosMarkers = zeros(N_cardinal_coor, N_markers);
    end
    
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

    for k = 1:N_markers
        trackpos = [model.markers.coordinates(:,k)' 1] * TransMatrix{model.markers.parent(k)}'; % Check for floating base?
        if isMX
            PosMarkers = {PosMarkers{:} trackpos(1:3)};
        else
            PosMarkers(:,k) = trackpos(1:3);
        end
    end
end