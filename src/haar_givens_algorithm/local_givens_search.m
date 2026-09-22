function [Q_best, f_best, info] = local_givens_search( ...
    Q0, obj, delta0, tol_step, max_sweeps, tol_improvement)
%LOCAL_GIVENS_SEARCH Coordinate search on SO(K) with Givens rotations.
%
%   At each sweep, all rotations Q*G_ij(+/-delta) are tested from the
%   current Q. The single best improving move is accepted. If no move
%   improves the objective, delta is halved.

Q_best = Q0;
f_best = obj(Q_best);
delta = delta0;
K = size(Q0, 1);
n_improvements = 0;
sweep = 0;

if nargin < 6 || isempty(tol_improvement)
    tol_improvement = 0;
end

for sweep = 1:max_sweeps
    if delta < tol_step
        break
    end

    best_trial_value = f_best;
    best_trial_Q = Q_best;

    for i = 1:K-1
        for j = i+1:K
            G_plus = local_givens_rotation(K, i, j, delta);
            Q_trial = Q_best * G_plus;
            f_trial = obj(Q_trial);

            if f_trial < best_trial_value - tol_improvement
                best_trial_value = f_trial;
                best_trial_Q = Q_trial;
            end

            G_minus = local_givens_rotation(K, i, j, -delta);
            Q_trial = Q_best * G_minus;
            f_trial = obj(Q_trial);

            if f_trial < best_trial_value - tol_improvement
                best_trial_value = f_trial;
                best_trial_Q = Q_trial;
            end
        end
    end

    if best_trial_value < f_best - tol_improvement
        Q_best = best_trial_Q;
        f_best = best_trial_value;
        n_improvements = n_improvements + 1;
    else
        delta = delta/2;
    end
end

info = struct();
info.sweeps = min(sweep, max_sweeps);
info.final_delta = delta;
info.n_improvements = n_improvements;
end
