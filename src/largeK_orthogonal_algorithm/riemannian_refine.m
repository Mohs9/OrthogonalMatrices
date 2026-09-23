function [Q, f, info] = riemannian_refine(Q, P, R, norm_func, opts)
%RIEMANNIAN_REFINE Projected subgradient refinement on O(K).
%
%   Minimizes a logarithmic version of the condition number while keeping Q
%   on the orthogonal group. Each iteration computes a tangent direction
%   Omega and updates with Q*expm(-alpha*Omega), which preserves
%   orthogonality without re-orthogonalization. R remains in the signature to
%   share the interface with the exact routines, although objective_smooth
%   uses P to reconstruct the inverse term.

if nargin < 4 || isempty(norm_func)
    norm_func = @norm_1;
end

% Allow riemannian_refine(Q, P, R, opts) without passing norm_func.
if nargin == 4 && isstruct(norm_func)
    opts = norm_func;
    norm_func = @norm_1;
end

if nargin < 5 || isempty(opts)
    opts = struct();
end

opts = fill_default_options_largeK(opts);
n_improvements = 0;

for iter = 1:opts.RIEMANN_MAX_ITERS
    [f, grad, ~, Omega] = objective_smooth(Q, P, norm_func);
    grad_norm = norm(grad, 'fro');

    % If the projected gradient norm is small, the iterate is near a
    % stationary point for the orthogonal geometry.
    if grad_norm < opts.RIEMANN_GRAD_TOL
        break
    end

    alpha = opts.ALPHA0;
    accepted = false;
    rel_improvement = 0;

    % Backtracking line search: shrink alpha until the smooth objective
    % strictly decreases.
    for iLine = 1:opts.MAX_LINESEARCH
        Qtrial = Q * expm(-alpha*Omega);
        ftrial = objective_smooth(Qtrial, P, norm_func);

        if ftrial < f
            rel_improvement = (f - ftrial)/max(1, abs(f));
            Q = Qtrial;
            f = ftrial;
            accepted = true;
            n_improvements = n_improvements + 1;
            break
        end

        alpha = alpha * opts.LINESEARCH_RHO;
    end

    % Stop if no acceptable step exists or the relative improvement is tiny.
    if ~accepted || rel_improvement < opts.RIEMANN_REL_TOL
        break
    end
end

% Recompute the final state to report value and gradient norm.
[f, grad] = objective_smooth(Q, P, norm_func);

info = struct();
info.iterations = iter;
info.n_improvements = n_improvements;
info.final_value = f;
info.final_grad_norm = norm(grad, 'fro');
end
