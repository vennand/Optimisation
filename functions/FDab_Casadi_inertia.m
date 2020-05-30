function [qdd] = FDab_Casadi_inertia(model, data, q, qd, tau, M, CoM, I)
    for i=1:length(data.segments)
        model.I{data.segments(i)} = mcI_Casadi(M{i},CoM{3*i-2:3*i}, diag(I{3*i-2:3*i}));
    end
    
    qdd = FDab_Casadi(model, q, qd, tau);
end