% ---------------------------------------------------------------------
% Global experiment for the large-K Haar + Riemannian + Givens algorithm.
%
% The script fixes one covariance matrix Sigma_e and runs
% optimize_Qstar_largeK many times with different random seeds. Each run
% stores the complete output struct returned by the optimizer.
% ---------------------------------------------------------------------

activate
tic

% Fixed reproducible covariance matrix.
rng(20, 'twister')
K = 30;
M = randn(K, K);
Sigma_e = M*M' + 0.25*eye(K);
P = chol(Sigma_e, 'lower');

% Number of independent searches for the same Sigma_e.
N = 50;

% Store outputs under a folder named after the matrix dimension.
plotDir = fullfile('plots', sprintf('K=%d', K));

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% Each cell stores one result struct returned by optimize_Qstar_largeK.
solutions = cell(N, 1);
kappa_values = zeros(N, 1);
improvement_values = zeros(N, 1);
objective_errors = zeros(N, 1);
orthogonality_errors = zeros(N, 1);
Q_best_values = cell(N, 1);

% Settings for each independent run.
opts = struct();
opts.SEED = 1;
opts.N_HAAR = 10000;
opts.N_ELITE = 30;
opts.RIEMANN_MAX_ITERS = 100;
opts.GIVENS_MAX_SWEEPS = 20;
opts.TOL_IMPROVEMENT = 1e-12;
seeds = opts.SEED + (0:N-1).';

% Objective norm used by the exact objective and by the Riemannian
% subgradient refinement.
norm_func = @norm_infinity;

parfor iRun = 1:N
    optsRun = opts;
    optsRun.SEED = seeds(iRun);

    result = optimize_Qstar_largeK(Sigma_e, norm_func, optsRun);

    solutions{iRun} = result;
    Q_best_values{iRun} = result.Q_star;
    kappa_values(iRun) = result.kappa_star;
    improvement_values(iRun) = result.improvement;
    objective_errors(iRun) = result.objective_error;
    orthogonality_errors(iRun) = result.orthogonality_error;

    fprintf('\nLarge-K Haar + Riemannian + Givens optimization\n')
    fprintf('K: %d\n', K)
    fprintf('Seed: %d\n', optsRun.SEED)
    fprintf('best Haar log(kappa_p): %.15g\n', result.best_haar_log_value)
    fprintf('best local kappa_p: %.15g\n', result.kappa_star)
    fprintf('kappa_p(P): %.15g\n', result.kappa_cholesky)
    fprintf('improvement: %.6f\n', result.improvement)
    fprintf('orthogonality error ||Q''Q-I||_F: %.3e\n', result.orthogonality_error)
    fprintf('det(Q_star): %.15g\n', result.det_Qstar)
    fprintf('cond(P*Q_star,p): %.15g\n', result.kappa_direct)
    fprintf('|cond(P*Q_star,p)-kappa_star|: %.3e\n', result.objective_error)

end

% Select the best solution among all independent runs.
[kappa_global, idx_best] = min(kappa_values);
result_global = solutions{idx_best};
Q_best_global = result_global.Q_star;
elapsed_time = toc;

fprintf('\nGlobal large-K experiment completed.\n')
fprintf('N: %d\n', N)
fprintf('K: %d\n', K)
fprintf('Total elapsed time: %.3f seconds\n', elapsed_time)
fprintf('Best run: %d\n', idx_best)
fprintf('kappa_p(P): %.12g\n', result_global.kappa_P)
fprintf('best kappa_p(P*Q): %.12g\n', kappa_global)
fprintf('best improvement: %.6f\n', result_global.improvement)
fprintf('best objective error: %.3e\n', result_global.objective_error)
fprintf('best orthogonality error: %.3e\n', result_global.orthogonality_error)

% Save the complete set of solutions and the inputs used in the experiment.
save(fullfile(plotDir, sprintf('largeK_orthogonal_K%d_solutions.mat', K)), ...
    'solutions', 'Q_best_values', 'kappa_values', 'improvement_values', ...
    'objective_errors', 'orthogonality_errors', ...
    'result_global', 'Q_best_global', ...
    'idx_best', 'P', 'Sigma_e', 'opts', 'seeds', 'norm_func', ...
    'N', 'elapsed_time')
