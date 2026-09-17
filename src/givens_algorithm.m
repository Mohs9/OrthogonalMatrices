function algorithm_results = givens_algorithm(P,  nRandom, nRefine)

K = size(P, 1);

if size(P, 2) ~= K
    error('P must be a square matrix.');
end

if rcond(P) < eps
    error('P must be nonsingular to compute its condition number.');
end


nTheta = K*(K-1)/2;
thetaLower = -pi*ones(1, nTheta);
thetaUpper =  pi*ones(1, nTheta);

objective = @(theta) condition_number( ...
    P*construct_Givens_matrices(theta, K), @norm_infinity);

%% Coarse exploration of the Givens angle space

% For K=2 this is a full one-dimensional grid over [-pi, pi].
% For larger K the Givens space has dimension K*(K-1)/2, so we use a
% random exploration and then refine the best candidates locally.

nGrid1D = 100000;

if nTheta == 1
    thetaCandidates = linspace(thetaLower, thetaUpper, nGrid1D).';
else
    U = lhsdesign(nRandom, nTheta);
    thetaCandidates = thetaLower + ...
        (thetaUpper-thetaLower).*U;
end

kappaCandidates = zeros(size(thetaCandidates, 1), 1);

for iCandidate = 1:size(thetaCandidates, 1)
    kappaCandidates(iCandidate) = objective(thetaCandidates(iCandidate, :));
end

[kappa_best, idxBest] = min(kappaCandidates);
theta_best = thetaCandidates(idxBest, :);
Q_best = construct_Givens_matrices(theta_best, K);

%% Local refinement from the best explored points

[~, idxSorted] = sort(kappaCandidates);
nRefine = min(nRefine, numel(idxSorted));

hasFmincon = exist('fmincon', 'file') == 2;

if hasFmincon
    options = optimoptions('fmincon', ...
        'Display', 'off', ...
        'Algorithm', 'active-set', ...
        'MaxIterations', 5000, ...
        'OptimalityTolerance', 1e-10, ...
        'StepTolerance', 1e-12);
else
    options = optimset( ...
        'Display', 'off', ...
        'MaxIter', 5000, ...
        'TolFun', 1e-10, ...
        'TolX', 1e-12);
end

kappas_nRefine = zeros(nRefine,1);
Q_nRefine =  cell(nRefine, 1);

for iStart = 1:nRefine
    theta0 = thetaCandidates(idxSorted(iStart), :);

    if hasFmincon
        [thetaCandidate, kappaCandidate] = fmincon( ...
            objective, theta0, ...
            [], [], [], [], ...
            thetaLower, thetaUpper, [], options);
    else
        unconstrainedObjective = @(theta) objective(wrap_to_pi_local(theta));
        [thetaCandidate, kappaCandidate] = fminsearch( ...
            unconstrainedObjective, theta0, options);
        thetaCandidate = wrap_to_pi_local(thetaCandidate);
    end

    kappas_nRefine(iStart) =  kappaCandidate;
    Q_nRefine{iStart} = construct_Givens_matrices(thetaCandidate,K);

    if kappaCandidate < kappa_best
        kappa_best = kappaCandidate;
        theta_best = thetaCandidate;
        Q_best = construct_Givens_matrices(theta_best, K);
    end
end

algorithm_results = struct();

algorithm_results.Q_best = Q_best;
algorithm_results.kappa_best = kappa_best;
algorithm_results.theta_best = theta_best;
algorithm_results.kappaCandidates = kappaCandidates;
algorithm_results.thetaCandidates = thetaCandidates;
algorithm_results.kappas_nRefine = kappas_nRefine;
algorithm_results.Q_nRefine = Q_nRefine;
end