function H = random_tangent_direction(K)
%RANDOM_TANGENT_DIRECTION Generate a normalized skew-symmetric direction.

H = zeros(K);

while norm(H, 'fro') == 0
    A = randn(K, K);
    H = A - A.';
end

H = H/norm(H, 'fro');
end