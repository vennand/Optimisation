function [qdd] = FDab_Casadi_gravity_inertia(model, q, qd, tau, G, M, CoM, I)
    pelvis = 6;
    thorax = 9;
    right_thigh = 33;
    left_thigh = 39;
    
%     model.I{pelvis} = mcI_Casadi(M(1),CoM(1,:),squeeze(I(1,:,:)));
%     model.I{thorax} = mcI_Casadi(M(2),CoM(2,:),squeeze(I(2,:,:)));
%     model.I{right_thigh} = mcI_Casadi(M(3),CoM(3,:),squeeze(I(3,:,:)));
%     model.I{left_thigh} = mcI_Casadi(M(4),CoM(4,:),squeeze(I(4,:,:)));
    
    model.I{pelvis} = mcI_Casadi(M(1),CoM(1,:),reshape(I(1,:),3,3)');
    model.I{thorax} = mcI_Casadi(M(2),CoM(2,:),reshape(I(2,:),3,3)');
    model.I{right_thigh} = mcI_Casadi(M(3),CoM(3,:),reshape(I(3,:),3,3)');
    model.I{left_thigh} = mcI_Casadi(M(4),CoM(4,:),reshape(I(4,:),3,3)');
    
    model.gravity = G;
    
    qdd = FDab_Casadi(model, q, qd, tau);
end