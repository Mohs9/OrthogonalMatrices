function Q = haar_SO(K)
%HAAR_SO Generate a Haar-distributed random matrix in SO(K).

Z = randn(K, K);
[Q, R] = qr(Z);

% Fix the sign ambiguity in the QR factorization.
d = sign(diag(R));
d(d == 0) = 1;
Q = Q * diag(d);

% Enforce determinant +1 by flipping one column if needed.
if det(Q) < 0
    Q(:,1) = -Q(:,1);
end
end
