function algorithm_results = givens_algorithm(P,  nRandom, nRefine, norma, opts)
%GIVENS_ALGORITHM Search for an orthogonal Givens factor that lowers kappa.
%
% Optional opts fields:
%   tracePaths      true saves theta/kappa values visited by local refiners.
%   initialization  'default' keeps the original K=2 grid behavior;
%                   'random' uses random starts for K=2.
%   nGrid1D         number of K=2 grid points for default initialization.
%   thetaCandidates user-supplied candidate starts, one row per start.
%   diagnosticTol   relative tolerance for active row diagnostics.
%   displaySummary  true displays [u, v] and active-row counts at the end.

if nargin < 5
    opts = struct();
end

if ~isfield(opts, 'tracePaths')
    opts.tracePaths = false;
end

if ~isfield(opts, 'initialization')
    opts.initialization = 'default';
end

if ~isfield(opts, 'nGrid1D')
    opts.nGrid1D = 100000;
end

if ~isfield(opts, 'thetaCandidates')
    opts.thetaCandidates = [];
end

if ~isfield(opts, 'diagnosticTol')
    opts.diagnosticTol = 1e-8;
end

if ~isfield(opts, 'displaySummary')
    opts.displaySummary = true;
end

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
    P*construct_Givens_matrices(theta, K), norma);

%% Coarse exploration of the Givens angle space

% For K=2 the default is a full one-dimensional grid over [-pi, pi].
% For larger K the Givens space has dimension K*(K-1)/2, so we use a
% random exploration and then refine the best candidates locally.

if ~isempty(opts.thetaCandidates)
    thetaCandidates = opts.thetaCandidates;
    if size(thetaCandidates, 2) ~= nTheta
        error('opts.thetaCandidates must have K*(K-1)/2 columns.');
    end
elseif nTheta == 1 && strcmpi(opts.initialization, 'random')
    thetaCandidates = thetaLower + ...
        (thetaUpper-thetaLower).*rand(nRandom, nTheta);
elseif nTheta == 1
    thetaCandidates = linspace(thetaLower, thetaUpper, opts.nGrid1D).';
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
thetaInitial_nRefine = zeros(nRefine, nTheta);
emptyResult = empty_solution_diagnostics();
local_results = repmat(emptyResult, nRefine, 1);

if opts.tracePaths
    thetaPaths = cell(nRefine, 1);
    kappaPaths = cell(nRefine, 1);
else
    thetaPaths = {};
    kappaPaths = {};
end

for iStart = 1:nRefine
    theta0 = thetaCandidates(idxSorted(iStart), :);
    thetaInitial_nRefine(iStart, :) = theta0;

    thetaPathCurrent = [];
    kappaPathCurrent = [];

    if hasFmincon
        if opts.tracePaths
            options_i = optimoptions(options, 'OutputFcn', @record_path);
        else
            options_i = options;
        end

        [thetaCandidate, kappaCandidate] = fmincon( ...
            objective, theta0, ...
            [], [], [], [], ...
            thetaLower, thetaUpper, [], options_i);
    else
        if opts.tracePaths
            options_i = optimset(options, 'OutputFcn', @record_path);
        else
            options_i = options;
        end

        unconstrainedObjective = @(theta) objective(wrap_to_pi_local(theta));
        [thetaCandidate, kappaCandidate] = fminsearch( ...
            unconstrainedObjective, theta0, options_i);
        thetaCandidate = wrap_to_pi_local(thetaCandidate);
    end

    if opts.tracePaths
        thetaPaths{iStart} = thetaPathCurrent;
        kappaPaths{iStart} = kappaPathCurrent;
    end

    kappas_nRefine(iStart) =  kappaCandidate;
    QCandidate = construct_Givens_matrices(thetaCandidate,K);
    Q_nRefine{iStart} = QCandidate;
    local_results(iStart) = solution_diagnostics(QCandidate);

    if kappaCandidate < kappa_best
        kappa_best = kappaCandidate;
        theta_best = thetaCandidate;
        Q_best = QCandidate;
    end
end

best_result = solution_diagnostics(Q_best);

if nRefine > 0
    final_results = [best_result; local_results];
else
    final_results = best_result;
end

algorithm_results = struct();

algorithm_results.Q_best = Q_best;
algorithm_results.kappa_best = kappa_best;
algorithm_results.theta_best = theta_best;
algorithm_results.kappaCandidates = kappaCandidates;
algorithm_results.thetaCandidates = thetaCandidates;
algorithm_results.kappas_nRefine = kappas_nRefine;
algorithm_results.Q_nRefine = Q_nRefine;
algorithm_results.thetaInitial_nRefine = thetaInitial_nRefine;
algorithm_results.thetaPaths = thetaPaths;
algorithm_results.kappaPaths = kappaPaths;
algorithm_results.best_result = best_result;
algorithm_results.local_results = local_results;
algorithm_results.results = final_results;

if opts.displaySummary
    for jResult = 1:numel(final_results)
        fprintf('\nRetained solution %d\n', jResult)
        disp('[u, v] =')
        disp([final_results(jResult).u, final_results(jResult).v])
        fprintf('numel(active_u): %d\n', numel(final_results(jResult).active_u))
        fprintf('numel(active_v): %d\n', numel(final_results(jResult).active_v))
    end
end

    function stop = record_path(theta, optimValues, state)
        stop = false;

        if strcmp(state, 'init') || strcmp(state, 'iter')
            thetaRow = theta(:).';

            if ~hasFmincon
                thetaRow = wrap_to_pi_local(thetaRow);
            end

            thetaPathCurrent(end+1, :) = thetaRow;

            if isfield(optimValues, 'fval')
                kappaPathCurrent(end+1, 1) = optimValues.fval;
            else
                kappaPathCurrent(end+1, 1) = objective(thetaRow);
            end
        end
    end

    function result = solution_diagnostics(Qstar)
        e = ones(K, 1);
        u = abs(P*Qstar)*e;
        v = abs(Qstar'/P)*e;
        max_u = max(u);
        max_v = max(v);
        active_u = find(abs(u-max_u) <= ...
            opts.diagnosticTol*max(1, abs(max_u)));
        active_v = find(abs(v-max_v) <= ...
            opts.diagnosticTol*max(1, abs(max_v)));

        result = struct();
        result.Q = Qstar;
        result.kappa = max_u*max_v;
        result.u = u;
        result.v = v;
        result.max_u = max_u;
        result.max_v = max_v;
        result.active_u = active_u;
        result.active_v = active_v;
    end

    function result = empty_solution_diagnostics()
        result = struct( ...
            'Q', [], ...
            'kappa', [], ...
            'u', [], ...
            'v', [], ...
            'max_u', [], ...
            'max_v', [], ...
            'active_u', [], ...
            'active_v', []);
    end
end
