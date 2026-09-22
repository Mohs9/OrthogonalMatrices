function P_matrices = permutation_matrices(K)
%PERMUTATION_MATRICES Generate all K-by-K permutation matrices.
%   P_matrices is a K! by 1 cell array. Each cell contains one
%   permutation matrix.

if K ~= floor(K) || K < 1
    error('K must be a positive integer.');
end

permutation_indices = perms(1:K);
nPermutations = size(permutation_indices, 1);
P_matrices = cell(nPermutations, 1);

for iPermutation = 1:nPermutations
    P = zeros(K);

    for iRow = 1:K
        P(iRow, permutation_indices(iPermutation, iRow)) = 1;
    end

    P_matrices{iPermutation} = P;
end
end
