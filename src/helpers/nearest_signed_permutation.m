function [distance_min, S_best, permutation_best, signs_best] = nearest_signed_permutation(Q1, Q2)
%NEAREST_SIGNED_PERMUTATION Find the closest signed permutation alignment.
%   [distance_min, S_best, permutation_best, signs_best] =
%   nearest_signed_permutation(Q1, Q2) solves
%
%       min_S ||Q2 - Q1*S||_F,
%
%   where S is a signed permutation matrix.
%
%   The signs are chosen from the entries of Q1'*Q2. The remaining problem
%   is a linear assignment over permutations, solved exactly by the
%   Hungarian algorithm in O(K^3), without enumerating 2^K K! candidates.

if ~isequal(size(Q1), size(Q2))
    error('Q1 and Q2 must have the same size.');
end

[nRows, K] = size(Q1);

if nRows ~= K
    error('Q1 and Q2 must be square.');
end

alignment = Q1.'*Q2;
score = abs(alignment);
permutation_best = hungarian_max_assignment(score);

signs_best = ones(1, K);
S_best = zeros(K);

for row = 1:K
    col = permutation_best(row);
    signValue = sign(alignment(row, col));

    if signValue == 0
        signValue = 1;
    end

    signs_best(row) = signValue;
    S_best(row, col) = signValue;
end

distance_min = norm(Q2 - Q1*S_best, 'fro');
end


