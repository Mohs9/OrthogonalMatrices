% ---------------------------------------------------------------------
% Inspect eigenvalues of products Q_i*Q_j from Givens solutions.
%
% Optional variables before running:
%   K           : dimension used to load plots/K=<K>/givens_K<K>_solutions.mat
%   productType : 'direct' for Q_i*Q_j, or 'relative' for Q_i'*Q_j
%
% The relative product Q_i'*Q_j is useful because it measures the rotation
% that maps Q_i into Q_j. Its eigenvalues lie on the unit circle when the
% matrices are orthogonal.
% ---------------------------------------------------------------------

clearvars -except solutions K plotDir productType
activate;

% Use relative rotations by default. Set productType = 'direct' before
% running the script to inspect Q_i*Q_j instead.
if ~exist('productType', 'var')
    productType = 'relative';
end

% Load solutions from disk unless they already exist in the workspace.
if ~exist('solutions', 'var')
    if ~exist('K', 'var')
        K = 3;
    end

    plotDir = fullfile('plots', sprintf('K=%d', K));
    solutionsFile = fullfile(plotDir, sprintf('givens_K%d_solutions.mat', K));
    load(solutionsFile, 'solutions')
else
    validIndex = find(~cellfun(@isempty, solutions), 1);
    K = size(solutions{validIndex}.Q_best, 1);

    if ~exist('plotDir', 'var')
        plotDir = fullfile('plots', sprintf('K=%d', K));
    end
end

% Create the output directory if needed.
if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% Remove empty cells, then determine the number of pairwise products.
validSolutions = ~cellfun(@isempty, solutions);
solutions = solutions(validSolutions);
N = numel(solutions);
nPairs = N*(N-1)/2;

% Preallocate storage for product matrices and their spectral summaries.
Q_products = zeros(K, K, nPairs);
eigenvaluesProduct = complex(zeros(nPairs, K));
eigenvalueAngles = zeros(nPairs, K);
eigenvalueModuli = zeros(nPairs, K);

pair_i = zeros(nPairs, 1);
pair_j = zeros(nPairs, 1);
detProduct = zeros(nPairs, 1);
traceProduct = zeros(nPairs, 1);
orthogonalityError = zeros(nPairs, 1);
maxModulusError = zeros(nPairs, 1);

idxPair = 0;

% Loop over all unordered pairs of solutions.
for i = 1:N-1
    Q_i = solutions{i}.Q_best;

    for j = i+1:N
        Q_j = solutions{j}.Q_best;
        idxPair = idxPair + 1;

        % Choose between direct products and relative rotations.
        switch productType
            case 'direct'
                Q_product = Q_i*Q_j;
            case 'relative'
                Q_product = Q_i.'*Q_j;
            otherwise
                error('productType must be ''direct'' or ''relative''.');
        end

        % Store eigenvalues and basic diagnostics for the product matrix.
        eigValues = eig(Q_product).';

        pair_i(idxPair) = i;
        pair_j(idxPair) = j;
        Q_products(:,:,idxPair) = Q_product;
        eigenvaluesProduct(idxPair,:) = eigValues;
        eigenvalueAngles(idxPair,:) = sort(angle(eigValues));
        eigenvalueModuli(idxPair,:) = abs(eigValues);
        detProduct(idxPair) = det(Q_product);
        traceProduct(idxPair) = trace(Q_product);
        orthogonalityError(idxPair) = norm(Q_product.'*Q_product - eye(K), 'fro');
        maxModulusError(idxPair) = max(abs(abs(eigValues)-1));
    end
end

% Print a short diagnostic summary.
fprintf('\nQ product eigenvalue analysis\n')
fprintf('K: %d\n', K)
fprintf('Solutions analyzed: %d\n', N)
fprintf('Pairs analyzed: %d\n', nPairs)
fprintf('Product type: %s\n', productType)
fprintf('Product matrices stored in Q_products: %d\n', size(Q_products, 3))
fprintf('Max orthogonality error: %.3e\n', max(orthogonalityError))
fprintf('Max eigenvalue modulus error from unit circle: %.3e\n', ...
    max(maxModulusError))

summaryTable = table( ...
    pair_i, ...
    pair_j, ...
    detProduct, ...
    traceProduct, ...
    orthogonalityError, ...
    maxModulusError);

% Save the full results. Q_products(:,:,k) corresponds to pair_i(k), pair_j(k).
outputStem = sprintf('Q_product_eigenvalues_%s', productType);
matFile = fullfile(plotDir, sprintf('%s.mat', outputStem));
summaryFile = fullfile(plotDir, sprintf('%s_summary.csv', outputStem));

save(matFile, ...
    'productType', ...
    'K', ...
    'N', ...
    'Q_products', ...
    'pair_i', ...
    'pair_j', ...
    'eigenvaluesProduct', ...
    'eigenvalueAngles', ...
    'eigenvalueModuli', ...
    'summaryTable')
writetable(summaryTable, summaryFile)

allEigenvalues = eigenvaluesProduct(:);
allAngles = eigenvalueAngles(:);
allModuli = eigenvalueModuli(:);

thetaCircle = linspace(0, 2*pi, 1000);

% Plot the distribution of eigenvalue angles.
figure
histogram(allAngles, 40)
xlabel('eigenvalue angle, radians')
ylabel('count')
title(sprintf('Eigenvalue angle distribution, %s products', productType))
grid on
saveas(gcf, fullfile(plotDir, sprintf('%s_angle_histogram.png', outputStem)))
