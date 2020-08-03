function angle = angle_between_vectors(v1, v2)
% angle = 2 * atan(norm(v1*norm(v2) - norm(v1)*v2) / ...
%     norm(v1 * norm(v2) + norm(v1) * v2))*180/pi;
angle = atan2(norm(cross(v1, v2)), dot(v1, v2))*180/pi;
end