function [Q_best, f_best, info] = local_tangent_refinement(Q0, obj, opts)
%LOCAL_TANGENT_REFINEMENT Try random tangent directions on SO(K).
%
%   Each accepted tangent perturbation is followed by a new Givens local
%   search. The tangent step is Q*expm(+/-epsilon*H), where H is skew
%   symmetric and normalized in Frobenius norm.

Q_best = Q0;
f_best = obj(Q_best);
K = size(Q0, 1);
epsilon = opts.TANGENT_RADIUS;

if ~isfield(opts, 'TOL_IMPROVEMENT') || isempty(opts.TOL_IMPROVEMENT)
    opts.TOL_IMPROVEMENT = 0;
end

rounds = 0;
n_tangent_improvements = 0;
n_givens_reruns = 0;
givens_after_tangent = repmat(struct( ...
    'sweeps', 0, ...
    'final_delta', NaN, ...
    'n_improvements', 0), opts.MAX_TANGENT_ROUNDS, 1);

while rounds < opts.MAX_TANGENT_ROUNDS && epsilon >= opts.TOL_TANGENT
    rounds = rounds + 1;
    best_trial_value = f_best;
    best_trial_Q = Q_best;

    for iDirection = 1:opts.N_TANGENT
        H = random_tangent_direction(K);

        Q_trial = Q_best * expm(epsilon*H);
        f_trial = obj(Q_trial);

        if f_trial < best_trial_value - opts.TOL_IMPROVEMENT
            best_trial_value = f_trial;
            best_trial_Q = Q_trial;
        end

        Q_trial = Q_best * expm(-epsilon*H);
        f_trial = obj(Q_trial);

        if f_trial < best_trial_value - opts.TOL_IMPROVEMENT
            best_trial_value = f_trial;
            best_trial_Q = Q_trial;
        end
    end

    if best_trial_value < f_best - opts.TOL_IMPROVEMENT
        Q_best = best_trial_Q;
        f_best = best_trial_value;
        n_tangent_improvements = n_tangent_improvements + 1;

        [Q_best, f_best, givens_info] = local_givens_search( ...
            Q_best, obj, opts.DELTA0, opts.TOL_STEP, opts.MAX_SWEEPS, ...
            opts.TOL_IMPROVEMENT);

        n_givens_reruns = n_givens_reruns + 1;
        givens_after_tangent(n_givens_reruns) = givens_info;
    else
        epsilon = epsilon/2;
    end
end

info = struct();
info.rounds = rounds;
info.final_epsilon = epsilon;
info.n_tangent_improvements = n_tangent_improvements;
info.n_givens_reruns = n_givens_reruns;
info.givens_after_tangent = givens_after_tangent(1:n_givens_reruns);
end
