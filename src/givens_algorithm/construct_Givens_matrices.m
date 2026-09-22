%% ================================================================
% Construct orthogonal matrix from Givens rotations
% ================================================================

function Q = construct_Givens_matrices(theta,K)
% construct_Givens_matrices constructs an orthogonal KxK matrix
% from Givens rotation angles.
%
% theta must have K*(K-1)/2 elements.

nTheta = K*(K-1)/2;

if length(theta) ~= nTheta
    error('theta must contain K*(K-1)/2 angles.');
end

Q = eye(K);

k = 1;

for i = 1:K-1
    for j = i+1:K

        c = cos(theta(k));
        s = sin(theta(k));

        G = eye(K);

        G(i,i) = c;
        G(i,j) = -s;
        G(j,i) = s;
        G(j,j) = c;

        Q = Q*G;

        k = k + 1;

    end
end

end