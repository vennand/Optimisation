function [data] = saveInitialValues(model, data)
N_segment = data.nSegment;
N_cardinal_coor = data.nCardinalCoor;

mass = zeros(N_segment,1);
CoM = zeros(N_segment,N_cardinal_coor);
I = zeros(N_segment,N_cardinal_coor);

segments = data.segments;

for i=1:N_segment
    [mass(i),CoM(i,:),tempI] = mcI(model.I{segments(i)});
    I(i,:) = diag(tempI);
end

data.initialMass = mass;
data.initialCoM = CoM;
data.initialInertia = I;
end