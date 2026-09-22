function g = norm_infinity(A)
%NORM_INFINITY Compute the induced infinity norm of a matrix.
%   For a matrix A, this norm is the maximum absolute row sum.

% Sum absolute values across each row.
rowSumA = sum(abs(A),2);

% Select the largest row sum.
g = max(rowSumA);


end
