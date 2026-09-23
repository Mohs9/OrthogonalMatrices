% Reproducible K=25 example for optimize_Qstar_largeK.
activate

rng(25, 'twister')

K = 20;
M = randn(K, K);
Sigma_e = M*M' + 0.25*eye(K);

opts = struct();
opts.SEED = 25;
opts.N_HAAR = 10000;
opts.N_ELITE = 30;
opts.RIEMANN_MAX_ITERS = 100;
opts.GIVENS_MAX_SWEEPS = 2;

result = optimize_Qstar_largeK(Sigma_e, @norm_infinity, opts);

fprintf('\nLarge-K Haar + Riemannian + Givens optimization\n')
fprintf('K: %d\n', K)

fprintf('best Haar log(kappa_p): %.15g\n', result.best_haar_log_value)
fprintf('best local kappa_p: %.15g\n', result.kappa_star)
fprintf('kappa_p(P): %.15g\n', result.kappa_cholesky)
fprintf('improvement: %.6f\n', result.improvement)
fprintf('orthogonality error ||Q''Q-I||_F: %.3e\n', result.orthogonality_error)
fprintf('det(Q_star): %.15g\n', result.det_Qstar)
fprintf('cond(P*Q_star,p): %.15g\n', result.kappa_direct)
fprintf('|cond(P*Q_star,p)-kappa_star|: %.3e\n', result.objective_error)


plotDir = fullfile('plots', sprintf('K=%d', K));

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

save(fullfile(plotDir, 'largeK_Qstar_result.mat'), 'result', 'Sigma_e', 'opts')
