run('startup.m')

% matrix = [0.9978916035 0.0226102511 0.0608368653 -0.0018212712;-0.0082063648 0.9737889385 -0.2273054305 0.0138886119;-0.0643816993 0.2263269310 0.9719213533 0.1596893073;0.0000000000 0.0000000000 0.0000000000 1.0000000000];
matrix = [matrice_rotation(rand(3,1),'xyz'), rand(3,1); 0, 0 , 0, 1];


tic
for i=1:10000
%     matrix = [matrice_rotation(rand(3,1),'xyz'), rand(3,1); 0, 0 , 0, 1];
    pluho(invR(matrix));
end
toc

tic
for i=1:10000
%     matrix = [matrice_rotation(rand(3,1),'xyz'), rand(3,1); 0, 0 , 0, 1];
    inv(pluho(matrix));
end
toc

pluho(invR(matrix))-inv(pluho(matrix))