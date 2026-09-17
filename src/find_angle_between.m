function angle = find_angle_between(vec1,vec2)

ab = dot(vec1,vec2);
a = norm(vec1);
b = norm(vec2);
angle = acos(ab/(a*b));

end