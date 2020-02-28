function [data] = saveInitialValues(model, data)
% total_mass = 0;
% for i=1:model.NB
%     [mass,~,~] = mcI(model.I{i});
%     
%     total_mass = total_mass + mass;
% end
% data.totalMass = total_mass;

if data.optimiseInertia
    N_segment = data.nSegment;
    N_cardinal_coor = data.nCardinalCoor;
    
    mass = zeros(N_segment,1);
    CoM = zeros(N_segment,N_cardinal_coor);
    I = zeros(N_segment,N_cardinal_coor,N_cardinal_coor);
    
    pelvis = 6;
    thorax = 9;
    right_thigh = 33;
    left_thigh = 39;
    
    segments = [pelvis, thorax, right_thigh, left_thigh];
    
    for i=1:N_segment
        [mass(i),CoM(i,:),I(i,:,:)] = mcI(model.I{segments(i)});
    end
    
    data.initialMass = mass;
    data.initialCoM = CoM;
    data.initialInertia = I;
end
end