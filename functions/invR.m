% Inverse 4x4 rotational-translation matrix
function out = invR(in)
R = in(1:3,1:3)';
out = [R, -R * in(1:3, 4); 0 0 0 1];
end