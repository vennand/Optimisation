function [qdd] = FDab_Casadi_inertia(model, q, qd, tau, I)
    model.I
    
    qdd = FDab_Casadi(model, q, qd, tau);
end