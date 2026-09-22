% ---------------------------------------------------------------------
% Global Haar + Givens experiment.
%
% The script fixes one covariance matrix Sigma_e and runs the
% Haar + Givens optimizer many times with different random seeds. Each run
% stores the complete output struct returned by haar_givens_algorithm.
% ---------------------------------------------------------------------

activate
tic

% Example covariance matrix Sigma_e = P*P'. It must be symmetric positive
% definite because haar_givens_algorithm starts with a Cholesky factor.
P = tril(randn(5,5));
Sigma_e = P*P.';
K = size(Sigma_e, 1);

% Number of independent searches for the same Sigma_e.
N = 50;

% Store outputs under a folder named after the matrix dimension.
plotDir = fullfile('plots', sprintf('K=%d', K));

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% Each cell stores one result struct returned by haar_givens_algorithm.
solutions = cell(N, 1);
kappa_values = zeros(N, 1);
improvement_values = zeros(N, 1);
Q_best_values = cell(N, 1);

% Settings for each independent run. These values favor deeper refinement
% per run instead of many shallow restarts.
opts = struct();
opts.SEED = 1;
opts.N_HAAR = 10000;
opts.N_ELITE = 30;
opts.MAX_SWEEPS = 1000;
opts.N_TANGENT = 100;
opts.MAX_TANGENT_ROUNDS = 15;
opts.TANGENT_RADIUS = 0.05;
opts.TOL_IMPROVEMENT = 1e-12;
seeds = opts.SEED + (0:N-1).';

parfor iRun = 1:N
    optsRun = opts;
    optsRun.SEED = seeds(iRun);

    result = haar_givens_algorithm(Sigma_e,@norm_infinity, optsRun);
    solutions{iRun} = result;
    Q_best_values{iRun} = result.Q_star;
    kappa_values(iRun) = result.kappa_star;
    improvement_values(iRun) = result.improvement;

    fprintf('\nHaar + Givens run %d/%d completed.\n', iRun, N)
    fprintf('seed: %d\n', optsRun.SEED)
    fprintf('best kappa_p(P*Q): %.12g\n', result.kappa_star)
    fprintf('improvement: %.6f\n', result.improvement)
    fprintf('orthogonality error ||Q''Q-I||_F: %.3e\n', result.orthogonality_error)
    disp('Q_best =')
    disp(result.Q_star)
end

% Select the best solution among all independent runs.
[kappa_global, idx_best] = min(kappa_values);
result_global = solutions{idx_best};
Q_best_global = result_global.Q_star;
elapsed_time = toc;

fprintf('\nGlobal Haar + Givens experiment completed.\n')
fprintf('N: %d\n', N)
fprintf('K: %d\n', K)
fprintf('Total elapsed time: %.3f seconds\n', elapsed_time)
fprintf('Best run: %d\n', idx_best)
fprintf('kappa_p(P): %.12g\n', result_global.kappa_P)
fprintf('best kappa_p(P*Q): %.12g\n', kappa_global)
fprintf('best improvement: %.6f\n', result_global.improvement)
disp('Q_best global =')
disp(Q_best_global)

% Save the complete set of solutions and the inputs used in the experiment.
save(fullfile(plotDir, sprintf('haar_givens_K%d_solutions.mat', K)), ...
    'solutions', 'Q_best_values', 'kappa_values', 'improvement_values', ...
    'result_global', 'Q_best_global', ...
    'idx_best', 'P', 'Sigma_e', 'opts', 'seeds', 'N', 'elapsed_time')
