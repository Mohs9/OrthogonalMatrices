function Q = orthogonal_matrix_generator(m, n)
%ORTHOGONAL_MATRIX_GENERATOR Generate a random matrix with orthonormal columns.
%   Q is obtained from a thin QR factorization of a Gaussian random matrix.
%   If n is omitted, the function returns an m-by-m orthogonal matrix.

% Use a square matrix by default.
if nargin < 2, n = m; end

% Compute the thin QR decomposition.
[Q, R] = qr(randn(m, n), 0);

% Adjust the signs to obtain the standard Haar-distributed QR convention.
d = diag(R);
ph = d ./ abs(d);
ph(isnan(ph)) = 1;
Q = Q * diag(ph);
end
