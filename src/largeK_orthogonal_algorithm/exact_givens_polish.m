function [Q, f, info] = exact_givens_polish(Q, P, R, norm_func, opts)
%EXACT_GIVENS_POLISH Final derivative-free polishing with Givens rotations.
%
%   Sweeps over the (i,j) planes of O(K) and, for each plane, searches a
%   small grid for the Givens angle that reduces the exact condition number
%   kappa(P*Q) = ||P*Q|| * ||Q'*R||. Since this phase uses objective_exact,
%   it validates the solution with the original objective after the smooth
%   Riemannian refinement.

if nargin < 4 || isempty(norm_func)
    norm_func = @norm_1;
end

% Allow exact_givens_polish(Q, P, R, opts) without passing norm_func.
if nargin == 4 && isstruct(norm_func)
    opts = norm_func;
    norm_func = @norm_1;
end

if nargin < 5 || isempty(opts)
    opts = struct();
end

opts = fill_default_options_largeK(opts);

K = size(Q, 1);
[f, ~] = objective_exact(Q, P, R, norm_func);
n_improvements = 0;
sweep_values = zeros(opts.GIVENS_MAX_SWEEPS, 1);

% Enumerate all Givens planes with i < j.
[plane_i, plane_j] = find(triu(ones(K), 1));
nPlanes = numel(plane_i);

for sweep = 1:opts.GIVENS_MAX_SWEEPS
    f_before = f;
    % Randomizing the order reduces dependence on a fixed plane ordering.
    order = randperm(nPlanes);

    for idx = order
        i = plane_i(idx);
        j = plane_j(idx);

        [theta_best, f_trial] = best_givens_angle(Q, P, R, i, j, norm_func, opts);

        % Accept only improvements larger than the numerical tolerance.
        if f_trial < f - opts.TOL_IMPROVEMENT
            Q = Q * local_givens_rotation(K, i, j, theta_best);
            f = f_trial;
            n_improvements = n_improvements + 1;
        end
    end

    sweep_values(sweep) = f;
    rel_improvement = (f_before - f)/max(1, abs(f_before));

    % Stop when a full sweep no longer gives a meaningful improvement.
    if rel_improvement < opts.TOL_GIVENS
        break
    end
end

% Store diagnostics useful for comparing elites and debugging convergence.
info = struct();
info.sweeps = sweep;
info.n_improvements = n_improvements;
info.sweep_values = sweep_values(1:sweep);
end
