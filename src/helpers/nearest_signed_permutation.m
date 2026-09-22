function [distance_min, S_best, permutation_best, signs_best] = nearest_signed_permutation(Q1, Q2)
%NEAREST_SIGNED_PERMUTATION Find the closest signed permutation alignment.
%   [distance_min, S_best, permutation_best, signs_best] =
%   nearest_signed_permutation(Q1, Q2) solves
%
%       min_S ||Q2 - Q1*S||_F,
%
%   where S is a signed permutation matrix.

if ~isequal(size(Q1), size(Q2))
    error('Q1 and Q2 must have the same size.');
end

[~, K] = size(Q1);

if K > 20
    error('This exact dynamic program uses 2^K states; K must be <= 20.');
end

alignment = Q1.'*Q2;
score = abs(alignment);
nMasks = 2^K;

bestScore = -Inf(nMasks, 1);
previousMask = zeros(nMasks, 1);
chosenColumn = zeros(nMasks, 1);
bestScore(1) = 0;

for mask = 0:nMasks-1
    row = sum(bitget(mask, 1:K)) + 1;

    if row > K || isinf(bestScore(mask+1))
        continue
    end

    for col = 1:K
        colBit = bitshift(1, col-1);

        if bitand(mask, colBit) ~= 0
            continue
        end

        newMask = bitor(mask, colBit);
        candidateScore = bestScore(mask+1) + score(row, col);

        if candidateScore > bestScore(newMask+1)
            bestScore(newMask+1) = candidateScore;
            previousMask(newMask+1) = mask;
            chosenColumn(newMask+1) = col;
        end
    end
end

permutation_best = zeros(1, K);
mask = nMasks-1;

for row = K:-1:1
    col = chosenColumn(mask+1);
    permutation_best(row) = col;
    mask = previousMask(mask+1);
end

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
