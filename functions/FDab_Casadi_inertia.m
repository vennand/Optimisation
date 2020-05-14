function [qdd] = FDab_Casadi_inertia(model, q, qd, tau, M, CoM, I)
    pelvis = 6;
    thorax = 9;
    right_thigh = 33;
    left_thigh = 39;
    
%     model.I{pelvis} = mcI_Casadi(M(1),CoM(1,:),squeeze(I(1,:,:)));
%     model.I{thorax} = mcI_Casadi(M(2),CoM(2,:),squeeze(I(2,:,:)));
%     model.I{right_thigh} = mcI_Casadi(M(3),CoM(3,:),squeeze(I(3,:,:)));
%     model.I{left_thigh} = mcI_Casadi(M(4),CoM(4,:),squeeze(I(4,:,:)));
    
    model.I{pelvis} = mcI_Casadi(M(1),CoM(1,:),diag(I(1,:)));
    model.I{thorax} = mcI_Casadi(M(2),CoM(2,:),diag(I(2,:)));
    model.I{right_thigh} = mcI_Casadi(M(3),CoM(3,:),diag(I(3,:)));
    model.I{left_thigh} = mcI_Casadi(M(4),CoM(4,:),diag(I(4,:)));
    
    qdd = FDab_Casadi(model, q, qd, tau);
end