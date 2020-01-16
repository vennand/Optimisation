% Calcule la matrice Rototrans passant des angles de cardans du bassin à [1000;0100;0010;0001]
% Pour Do_822_contact_2_MOD200.00_GenderF_DoCig
function [Roto] = refential_matrix()

angleX = 0.0480;
angleY = -0.0657;
angleZ = 1.5720;

RotX = [1, 0, 0;...
        0, cos(angleX), -sin(angleX);...
        0, sin(angleX), cos(angleX)];


RotY = [cos(angleY), 0, sin(angleY);...
        0, 1, 0;...
        -sin(angleY), 0, cos(angleY)];

    
RotZ = [cos(angleZ), -sin(angleZ), 0;...
        sin(angleZ), cos(angleZ), 0;...
        0, 0, 1];

Roto = RotX  * RotY * RotZ;
end



