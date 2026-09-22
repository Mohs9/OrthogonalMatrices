function result = optimize_Qstar_largeK(Sigma_e, norm_func, opts)
%OPTIMIZE_QSTAR_LARGEK Haar + Riemannian + Givens search for large K.
%
%   This is the large-K analogue of haar_givens_algorithm:
%     1. Haar exploration on O(K).
%     2. Keep the best elites.
%     3. Refine each elite with a smooth Riemannian descent.
%     4. Polish with exact 1-norm Givens rotations.

if nargin < 2 || isempty(norm_func)
    norm_func = @norm_1;
end

if nargin == 2 && isstruct(norm_func)
    opts = struct();
    opts = norm_func;
    norm_func = @norm_1;
end

if nargin < 3 || isempty(opts)
    opts = struct();
end

opts = fill_default_options_largeK(opts);

if ~isempty(opts.SEED)
    rng(opts.SEED, 'twister');
end

if size(Sigma_e, 1) ~= size(Sigma_e, 2)
    error('Sigma_e must be square.');
end

Sigma_e = (Sigma_e + Sigma_e')/2;
K = size(Sigma_e, 1);
P = chol(Sigma_e, 'lower');
R = P \ eye(K);

% ---------------------------------------------------------------------
% 1. Global Haar exploration on O(K).
% ---------------------------------------------------------------------
nCandidates = opts.N_HAAR + 1;
haar_values = zeros(nCandidates, 1);
Q_haar = zeros(K, K, nCandidates);

Q_haar(:,:,1) = eye(K);
[~, haar_values(1)] = objective_exact(Q_haar(:,:,1), P, R, norm_func);

for i = 2:nCandidates
    Q = orthogonal_matrix_generator(K);
    Q_haar(:,:,i) = Q;
    [~, haar_values(i)] = objective_exact(Q, P, R, norm_func);
end

% ---------------------------------------------------------------------
% 2. Elite selection.
% ---------------------------------------------------------------------
[haar_values_sorted, idx_sorted] = sort(haar_values, 'ascend');
nElite = min(opts.N_ELITE, nCandidates);
elite_indices = idx_sorted(1:nElite);
best_haar_log_value = haar_values_sorted(1);

% ---------------------------------------------------------------------
% 3. Smooth Riemannian refinement + exact Givens polishing.
% ---------------------------------------------------------------------
local_values = zeros(nElite, 1);
Q_local = zeros(K, K, nElite);
local_info = repmat(struct('riemannian', [], 'givens', []), nElite, 1);

for iElite = 1:nElite
    Q0 = Q_haar(:,:,elite_indices(iElite));

    [Q_riem, ~, riem_info] = riemannian_refine(Q0, P, R, norm_func, opts);
    [Q_best, kappa_best, givens_info] = exact_givens_polish( ...
        Q_riem, P, R, norm_func, opts);

    Q_local(:,:,iElite) = Q_best;
    local_values(iElite) = kappa_best;
    local_info(iElite).riemannian = riem_info;
    local_info(iElite).givens = givens_info;
end

% ---------------------------------------------------------------------
% 4. Final solution and checks.
% ---------------------------------------------------------------------
[kappa_star, idx_best_local] = min(local_values);
Q_star = Q_local(:,:,idx_best_local);
B0_inv_star = P * Q_star;

orthogonality_error = norm(Q_star'*Q_star - eye(K), 'fro');
det_Qstar = det(Q_star);
[kappa_direct, ~] = objective_exact(Q_star, P, R, norm_func);
objective_error = abs(kappa_direct - kappa_star);
[kappa_P, ~] = objective_exact(eye(K), P, R, norm_func);
improvement = (kappa_P - kappa_star)/kappa_P;

fprintf('\nLarge-K Haar + Riemannian + Givens optimization\n')
fprintf('K: %d\n', K)
fprintf('N_HAAR: %d\n', opts.N_HAAR)
fprintf('N_ELITE: %d\n', nElite)
fprintf('best Haar log(kappa_p): %.15g\n', best_haar_log_value)
fprintf('best local kappa_p: %.15g\n', kappa_star)
fprintf('kappa_p(P): %.15g\n', kappa_P)
fprintf('improvement: %.6f\n', improvement)
fprintf('orthogonality error ||Q''Q-I||_F: %.3e\n', orthogonality_error)
fprintf('det(Q_star): %.15g\n', det_Qstar)
fprintf('cond(P*Q_star,p): %.15g\n', kappa_direct)
fprintf('|cond(P*Q_star,p)-kappa_star|: %.3e\n', objective_error)

result = struct();
result.Q_star = Q_star;
result.kappa_star = kappa_star;
result.B0_inv_star = B0_inv_star;
result.best_haar_log_value = best_haar_log_value;
result.local_values = local_values;
result.Q_local = Q_local;
result.local_info = local_info;
result.P = P;
result.R = R;
result.kappa_P = kappa_P;
result.kappa_cholesky = kappa_P;
result.improvement = improvement;
result.orthogonality_error = orthogonality_error;
result.det_Qstar = det_Qstar;
result.det_Q_star = det_Qstar;
result.kappa_direct = kappa_direct;
result.objective_error = objective_error;
result.haar_values = haar_values;
result.elite_indices = elite_indices;
result.opts = opts;
result.norm_func = func2str(norm_func);
end
