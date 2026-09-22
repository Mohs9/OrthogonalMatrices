% Reproducible K=25 example for optimize_Qstar_largeK.
activate

rng(25, 'twister')

K = 25;
M = randn(K, K);
Sigma_e = M*M' + 0.25*eye(K);

opts = struct();
opts.SEED = 25;
opts.N_HAAR = 10000;
opts.N_ELITE = 10;
opts.RIEMANN_MAX_ITERS = 100;
opts.GIVENS_MAX_SWEEPS = 2;

result = optimize_Qstar_largeK(Sigma_e, @norm_infinity, opts);

fprintf('\nK = %d\n', K)
fprintf('kappa(P): %.15g\n', result.kappa_cholesky)
fprintf('kappa_star: %.15g\n', result.kappa_star)
fprintf('relative improvement: %.6f\n', result.improvement)
fprintf('N_HAAR: %d\n', result.opts.N_HAAR)
fprintf('N_ELITE: %d\n', result.opts.N_ELITE)
fprintf('orthogonality error: %.3e\n', result.orthogonality_error)
fprintf('objective error: %.3e\n', result.objective_error)

plotDir = fullfile('plots', sprintf('K=%d', K));

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

save(fullfile(plotDir, 'largeK_Qstar_result.mat'), 'result', 'Sigma_e', 'opts')
