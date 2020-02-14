function [qdd] = FDab_Casadi_gravity(model, q, qd, tau, G)
    model.gravity = G;
    
    qdd = FDab_Casadi(model, q, qd, tau);
end