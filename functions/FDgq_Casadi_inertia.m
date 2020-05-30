function [qdd] = FDgq_Casadi_inertia(model, data, q, qd, tau, M, CoM, I)
    for i=1:length(data.segments)
        model.I{data.segments(i)} = mcI_Casadi(M(i),CoM(i,:), diag(I(i,:)));
    end
    
    qdd = FDgq_Casadi(model, q, qd, tau);
end