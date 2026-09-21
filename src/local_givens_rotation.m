function G = local_givens_rotation(K, i, j, theta)
%LOCAL_GIVENS_ROTATION Construct a local Givens rotation in the (i,j) plane.

if ~(1 <= i && i < j && j <= K)
    error('Indices must satisfy 1 <= i < j <= K.');
end

c = cos(theta);
s = sin(theta);

G = eye(K);
G(i,i) = c;
G(i,j) = -s;
G(j,i) = s;
G(j,j) = c;
end
