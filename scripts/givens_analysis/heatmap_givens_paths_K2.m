% ---------------------------------------------------------------------
% Heat map of local Givens-refinement paths for K=2.
%
% For K=2 there is only one Givens angle. This script starts the local
% optimizer from many random initial angles and shows where the optimizer
% spends its iterations.
% ---------------------------------------------------------------------

activate
rng(10, 'twister')

K = 2;
nStarts = 1000;
nBinsTheta = 200;
plotDir = fullfile('plots', sprintf('K=%d', K));

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% Reuse an existing 2-by-2 P from the workspace if available. Otherwise,
% generate one lower-triangular test matrix.
if ~exist('P', 'var') || ~isequal(size(P), [K K])
    P = tril(randn(K,K));

    while rcond(P) < eps
        P = tril(randn(K,K));
    end
end

opts = struct();
opts.initialization = 'random';
opts.tracePaths = true;

results = givens_algorithm(P, nStarts, nStarts, @norm_infinity, opts);

allTheta = [];
allIter = [];

for iStart = 1:numel(results.thetaPaths)
    thetaPath = results.thetaPaths{iStart};

    if isempty(thetaPath)
        continue
    end

    nPath = size(thetaPath, 1);
    allTheta = [allTheta; thetaPath(:, 1)];
    allIter = [allIter; (1:nPath).'];
end

maxIter = max(allIter);
thetaEdges = linspace(-pi, pi, nBinsTheta + 1);
iterEdges = 0.5:1:(maxIter + 0.5);

counts = histcounts2(allTheta, allIter, thetaEdges, iterEdges);
thetaCenters = 0.5*(thetaEdges(1:end-1) + thetaEdges(2:end));
iterCenters = 1:maxIter;

figure
imagesc(thetaCenters, iterCenters, counts.')
axis xy
xlim([-pi pi])
xlabel('\theta')
ylabel('Iteracion')
title('Recorrido del algoritmo desde condiciones iniciales aleatorias')
cb = colorbar;
ylabel(cb, 'Visitas')
grid on
saveas(gcf, fullfile(plotDir, 'givens_K2_path_heatmap.png'))
savefig(gcf, fullfile(plotDir, 'givens_K2_path_heatmap.fig'))

figure
plot(results.thetaCandidates, results.kappaCandidates, '.', ...
    'MarkerSize', 8)
hold on
plot(results.theta_best, results.kappa_best, 'ro', ...
    'MarkerFaceColor', 'r')
xlabel('\theta inicial')
ylabel('\kappa_\infty(PQ)')
title('Condiciones iniciales y mejor solucion encontrada')
grid on
saveas(gcf, fullfile(plotDir, 'givens_K2_initial_conditions.png'))
savefig(gcf, fullfile(plotDir, 'givens_K2_initial_conditions.fig'))

fprintf('\nK=2 Givens path heat map completed.\n')
fprintf('Number of starts: %d\n', nStarts)
fprintf('Figures saved in: %s\n', plotDir)
fprintf('kappa_inf(P): %.12g\n', condition_number(P, @norm_infinity))
fprintf('best kappa_inf(P*Q): %.12g\n', results.kappa_best)
fprintf('theta_best degrees: %.12g\n', results.theta_best*(180/pi))
