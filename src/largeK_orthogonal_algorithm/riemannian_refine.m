function [Q, f, info] = riemannian_refine(Q, P, R, norm_func, opts)
%RIEMANNIAN_REFINE Projected subgradient descent on O(K).

if nargin < 4 || isempty(norm_func)
    norm_func = @norm_1;
end

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

    if grad_norm < opts.RIEMANN_GRAD_TOL
        break
    end

    alpha = opts.ALPHA0;
    accepted = false;
    rel_improvement = 0;

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

    if ~accepted || rel_improvement < opts.RIEMANN_REL_TOL
        break
    end
end

[f, grad] = objective_smooth(Q, P, norm_func);

info = struct();
info.iterations = iter;
info.n_improvements = n_improvements;
info.final_value = f;
info.final_grad_norm = norm(grad, 'fro');
end
