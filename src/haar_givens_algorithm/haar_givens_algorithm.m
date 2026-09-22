function result = haar_givens_algorithm(Sigma_e, norm_func, opts)
%HAAR_GIVENS_ALGORITHM Minimize cond_1(P*Q) over Q in SO(K).
%
%   result = haar_givens_algorithm(Sigma_e, opts) computes the lower
%   Cholesky factor P such that Sigma_e = P*P' and searches for
%
%       Q_star = argmin_{Q in SO(K)} cond_p(P*Q).
%
%   The algorithm is deliberately simple:
%     1. Haar exploration on SO(K).
%     2. Keep the best N_ELITE matrices.
%     3. Refine each elite with explicit Givens rotations.
%     4. Try random tangent perturbations and re-run Givens after
%        successful perturbations.

if nargin < 2 || isempty(norm_func)
    norm_func = @norm_1;
end

if nargin == 2 && isstruct(norm_func)
    opts = norm_func;
    norm_func = @norm_1;
end

if nargin < 3 || isempty(opts)
    opts = struct();
end

if size(Sigma_e, 1) ~= size(Sigma_e, 2)
    error('Sigma_e must be square.');
end

K = size(Sigma_e, 1);
Sigma_e = (Sigma_e + Sigma_e.')/2;
P = chol(Sigma_e, 'lower');
P_inv = P \ eye(K);

opts = fill_default_options(opts);

if ~isempty(opts.SEED)
    rng(opts.SEED, 'twister');
end

obj = @(Q) condition_number(P*Q, norm_func);

% ---------------------------------------------------------------------
% 1. Global Haar exploration on SO(K).
% ---------------------------------------------------------------------
haar_values = zeros(opts.N_HAAR, 1);
Q_haar = zeros(K, K, opts.N_HAAR);

for i = 1:opts.N_HAAR
    Q = haar_SO(K);
    Q_haar(:,:,i) = Q;
    haar_values(i) = obj(Q);
end

% ---------------------------------------------------------------------
% 3. Elite selection.
% ---------------------------------------------------------------------
[haar_values_sorted, idx_sorted] = sort(haar_values, 'ascend');
nElite = min(opts.N_ELITE, opts.N_HAAR);
elite_indices = idx_sorted(1:nElite);
best_haar_value = haar_values_sorted(1);

% ---------------------------------------------------------------------
% 4. Local refinement with explicit Givens rotations and tangent directions.
% ---------------------------------------------------------------------
local_values = zeros(nElite, 1);
Q_local = zeros(K, K, nElite);
local_info = repmat(struct( ...
    'initial_givens', [], ...
    'tangent', []), nElite, 1);

for iElite = 1:nElite
    Q0 = Q_haar(:,:,elite_indices(iElite));
    [Q_givens, f_givens, givens_info] = local_givens_search( ...
        Q0, obj, opts.DELTA0, opts.TOL_STEP, opts.MAX_SWEEPS, ...
        opts.TOL_IMPROVEMENT);

    [Q_best, f_best, tangent_info] = local_tangent_refinement( ...
        Q_givens, obj, opts);

    Q_local(:,:,iElite) = Q_best;
    local_values(iElite) = f_best;
    local_info(iElite).initial_givens = givens_info;
    local_info(iElite).tangent = tangent_info;
end

% ---------------------------------------------------------------------
% 5. Final solution and numerical checks.
% ---------------------------------------------------------------------
[kappa_star, idx_best_local] = min(local_values);
Q_star = Q_local(:,:,idx_best_local);
B0_inv_star = P * Q_star;

orthogonality_error = norm(Q_star.'*Q_star - eye(K), 'fro');
det_Q_star = det(Q_star);
cond_direct = condition_number(B0_inv_star, norm_func);
cond_error = abs(cond_direct - kappa_star);
kappa_P = condition_number(P,norm_func);
improvement = (kappa_P - kappa_star)/kappa_P;

fprintf('\nHaar + Givens SO(K) optimization\n')
fprintf('K: %d\n', K)
fprintf('N_HAAR: %d\n', opts.N_HAAR)
fprintf('N_ELITE: %d\n', nElite)
fprintf('N_TANGENT: %d\n', opts.N_TANGENT)
fprintf('TANGENT_RADIUS: %.3g\n', opts.TANGENT_RADIUS)
fprintf('MAX_TANGENT_ROUNDS: %d\n', opts.MAX_TANGENT_ROUNDS)
fprintf('best Haar kappa_p: %.15g\n', best_haar_value)
fprintf('best local kappa_p: %.15g\n', kappa_star)
fprintf('kappa_p(P): %.15g\n', kappa_P)
fprintf('improvement: %.6f\n', improvement)
fprintf('orthogonality error ||Q''Q-I||_F: %.3e\n', orthogonality_error)
fprintf('det(Q_star): %.15g\n', det_Q_star)
fprintf('cond(P*Q_star,p): %.15g\n', cond_direct)
fprintf('|cond(P*Q_star,p)-kappa_star|: %.3e\n', cond_error)

result = struct();
result.Q_star = Q_star;
result.kappa_star = kappa_star;
result.B0_inv_star = B0_inv_star;
result.best_haar_value = best_haar_value;
result.local_values = local_values;
result.Q_local = Q_local;
result.local_info = local_info;
result.P = P;
result.P_inv = P_inv;
result.kappa_P = kappa_P;
result.improvement = improvement;
result.orthogonality_error = orthogonality_error;
result.det_Q_star = det_Q_star;
result.det_Qstar = det_Q_star;
result.cond_direct = cond_direct;
result.cond_error = cond_error;
result.cond1_direct = cond_direct;
result.cond1_error = cond_error;
result.haar_values = haar_values;
result.elite_indices = elite_indices;
result.opts = opts;
end
