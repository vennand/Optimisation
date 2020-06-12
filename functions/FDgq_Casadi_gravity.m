function [qdd] = FDgq_Casadi_gravity(model, q, qd, tau, G)
    model.gravity = G;
    
    qdd = FDgq_Casadi(model, q, qd, tau);
end