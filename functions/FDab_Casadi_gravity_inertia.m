function [qdd] = FDab_Casadi_gravity_inertia(model, q, qd, tau, G, M, CoM, I)
    pelvis = 6;
    thorax = 9;
    right_thigh = 33;
    left_thigh = 39;
    
    model.I{pelvis} = mcI([M(1),CoM(1,:),I(1,:,:)]);
    model.I{thorax} = mcI([M(2),CoM(2,:),I(2,:,:)]);
    model.I{right_thigh} = mcI([M(3),CoM(3,:),I(3,:,:)]);
    model.I{left_thigh} = mcI([M(4),CoM(4,:),I(4,:,:)]);
    
    model.gravity = G;
    
    qdd = FDab_Casadi(model, q, qd, tau);
end