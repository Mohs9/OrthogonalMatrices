% Simple smoke test for optimize_Qstar_largeK.

activate

rng(11, 'twister')

K = 6;
M = randn(K, K);
Sigma_e = M*M' + 0.5*eye(K);

opts = struct();
opts.SEED = 11;
opts.N_HAAR = 30;
opts.N_ELITE = 3;
opts.RIEMANN_MAX_ITERS = 5;
opts.GIVENS_MAX_SWEEPS = 1;

result = optimize_Qstar_largeK(Sigma_e, @norm_infinity, opts);

assert(isfield(result, 'Q_star'))
assert(isfield(result, 'kappa_star'))
assert(norm(result.Q_star'*result.Q_star - eye(K), 'fro') < 1e-8)
[kappa_check, ~] = objective_exact(result.Q_star, result.P, result.R, @norm_infinity);
assert(abs(kappa_check - result.kappa_star) < 1e-8)

fprintf('\nSimple large-K optimizer test passed.\n')
