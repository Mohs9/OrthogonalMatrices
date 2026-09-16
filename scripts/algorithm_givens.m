% ---------------------------------------------------------------------
% Explore the space of orthogonal matrices using Givens rotations.
%
% The objective is
%
%   min_Q kappa_inf(P*Q),
%
% where Q is built from K*(K-1)/2 Givens angles.
% ---------------------------------------------------------------------
activate
rng(1, 'twister')
tic


%% Problem setup
P = tril(randn(2,2));


% If P already exists in the workspace, the script uses it. Otherwise it
% starts with a small reproducible example.
if ~exist('P', 'var')
    P = [5, 0,;
        1, 6];
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
    P*construct_Givens_matrices(theta, K), @norm_infinity);

%% Coarse exploration of the Givens angle space

% For K=2 this is a full one-dimensional grid over [-pi, pi].
% For larger K the Givens space has dimension K*(K-1)/2, so we use a
% random exploration and then refine the best candidates locally.
nGrid1D = 10000;
nRandom = 100000;
nRefine = 100;

if nTheta == 1
    thetaCandidates = linspace(thetaLower, thetaUpper, nGrid1D).';
else
    thetaCandidates = thetaLower + ...
        (thetaUpper-thetaLower).*rand(nRandom, nTheta);
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
        'Algorithm', 'sqp', ...
        'MaxIterations', 1000, ...
        'OptimalityTolerance', 1e-10, ...
        'StepTolerance', 1e-12);
else
    options = optimset( ...
        'Display', 'off', ...
        'MaxIter', 5000, ...
        'TolFun', 1e-10, ...
        'TolX', 1e-12);
end

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

    if kappaCandidate < kappa_best
        kappa_best = kappaCandidate;
        theta_best = thetaCandidate;
        Q_best = construct_Givens_matrices(theta_best, K);
    end
end

%% Diagnostics

kappa_identity = condition_number(P, @norm_infinity);
orthogonality_error = norm(Q_best.'*Q_best - eye(K), 'fro');
elapsed_time = toc;

fprintf('\nGivens search completed.\n')
fprintf('K: %d\n', K)
fprintf('Number of Givens angles: %d\n', nTheta)
fprintf('kappa_inf(P): %.12g\n', kappa_identity)
fprintf('best kappa_inf(P*Q): %.12g\n', kappa_best)
fprintf('orthogonality error ||Q''Q-I||_F: %.3e\n', orthogonality_error)
fprintf('elapsed time: %.3f seconds\n\n', elapsed_time)

disp('theta_best =')
disp(theta_best*(360/(2*pi)))

disp('Q_best =')
disp(Q_best)

if nTheta == 1
    figure
    plot(thetaCandidates, kappaCandidates, 'LineWidth', 1.2)
    hold on
    plot(theta_best, kappa_best, 'ro', 'MarkerFaceColor', 'r')
    xlabel('\theta')
    ylabel('\kappa_\infty(PQ)')
    title('Givens angle exploration')
    grid on
elseif nTheta == 3
    figure
    scatter3( ...
        thetaCandidates(:, 1), ...
        thetaCandidates(:, 2), ...
        thetaCandidates(:, 3), ...
        18, log10(kappaCandidates), 'filled')
    hold on
    scatter3(theta_best(1), theta_best(2), theta_best(3), ...
        90, log10(kappa_best), 'r', 'filled')
    xlabel('\theta_{12}')
    ylabel('\theta_{13}')
    zlabel('\theta_{23}')
    title('Givens angle exploration for K=3')
    cb = colorbar;
    ylabel(cb, 'log_{10}(\kappa_\infty(PQ))')
    grid on
    view(45, 25)
end
