% Group Q_star solutions up to signed permutations.
%
% Two matrices Q_i and Q_j are considered equivalent if there exists a
% signed permutation matrix S such that
%
%   Q_j ~= Q_i*S.
%
% The distance used below is
%
%   min_S ||Q_j - Q_i*S||_F.

activate

% Load the stored solutions for the selected dimension.
K = 5;
plotDir = fullfile('plots', sprintf('K=%d', K));
load(fullfile(plotDir, sprintf('haar_givens_K%d_solutions.mat', K)))

N = numel(solutions);

% Tolerance for deciding whether the signed-permutation distance is zero
% numerically.
tol = 1e-5;

% D_perm(i,j) stores the distance between Q_i and Q_j after the best signed
% permutation alignment.
D_perm = zeros(N);

% Compute the pairwise signed-permutation distance matrix.
for i = 1:N
    Q_i = solutions{i}.Q_star;

    for j = i+1:N
        Q_j = solutions{j}.Q_star;

        D_perm(i,j) = nearest_signed_permutation(Q_i, Q_j);
        D_perm(j,i) = D_perm(i,j);
    end
end

% same_perm(i,j) is true when Q_i and Q_j are equivalent up to signed
% permutation under the selected tolerance.
same_perm = D_perm <= tol;

% Build connected groups of equivalent solutions.
% If Q_1 is equivalent to Q_2 and Q_2 is equivalent to Q_3, then the three
% are assigned to the same group.
groups = zeros(N,1);
g = 0;

for i = 1:N
    if groups(i) ~= 0
        continue
    end

    g = g + 1;
    stack = i;
    groups(i) = g;

    while ~isempty(stack)
        current = stack(end);
        stack(end) = [];

        neighbors = find(same_perm(current,:));
        for nb = neighbors
            if groups(nb) == 0
                groups(nb) = g;
                stack(end+1) = nb; %#ok<SAGROW>
            end
        end
    end
end

% Count how many solutions fall in each equivalence group.
nGroups = max(groups);
group_sizes = accumarray(groups, 1);

% Count equivalent pairs, excluding the diagonal.
nPairs = N*(N-1)/2;
nSamePairs = nnz(triu(same_perm, 1));

% Print summary statistics.
fprintf('\nSigned permutation groups\n')
fprintf('K: %d\n', K)
fprintf('Tolerance: %.1e\n', tol)
fprintf('Total solutions: %d\n', N)
fprintf('Equivalent pairs: %d of %d\n', nSamePairs, nPairs)
fprintf('Number of groups: %d\n', nGroups)

group_table = table((1:nGroups).', group_sizes, ...
    'VariableNames', {'Group', 'Size'});

% Display the size of each equivalence group.
disp(group_table)
