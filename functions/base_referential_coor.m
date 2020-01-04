% From the coordinates of a point in a segment referential, finds the coordinates in the base referential

function PosMarkers = base_referential_coor(model, q)
    import casadi.*
    
    N = model.NB;

%     Xa = cell(N,1);
%     TransMatrix = cell(N,1);

    isMX = false;

    N_markers = size(model.markers.coordinates);
    if class(q) == "casadi.SX"
        PosMarkers = SX.sym('marker', N_markers(1) * N_markers(2), 1);
    elseif class(q) == "casadi.MX"
        PosMarkers = {};
        isMX = true;
    elseif class(q) == "casadi.DM"
        PosMarkers = DM(N_markers(1) * N_markers(2), 1);
    else
        PosMarkers = zeros(N_markers(1) * N_markers(2), 1);
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

    n = 0;
    for k = 1:N_markers(2)
        trackpos = [model.markers.coordinates(:,k)' 1] * TransMatrix{model.markers.parent(k)}'; % Check for floating base?
        for m = 1:N_markers(1)
            n = n + 1;
            if isMX
                PosMarkers = {PosMarkers{:} trackpos(m)};
            else
                PosMarkers(n) = trackpos(m);
            end
        end
    end
end