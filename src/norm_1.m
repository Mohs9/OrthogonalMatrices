function g = norm_1(A)
%NORM_1 Compute the induced one norm of a matrix.
%   For a matrix A, this norm is the maximum absolute row sum.

% Sum absolute values across each row.
colSumA = sum(abs(A),1);

% Select the largest row sum.
g = max(colSumA);


end
