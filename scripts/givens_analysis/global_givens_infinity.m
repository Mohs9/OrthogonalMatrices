% ---------------------------------------------------------------------
% Global Givens search for the infinity-norm condition number objective.
%
% The script fixes one matrix P and runs the Givens-angle optimizer from
% many Latin-hypercube initial designs. Each run stores the best Q found,
% so the resulting cell array can later be used to analyze repetitions,
% signed-permutation equivalence, eigenvalues, and angle patterns.
% ---------------------------------------------------------------------

activate
%rng(1, 'twister')
tic

% Number of independent searches to run for the same matrix P.
N = 2;

% Generate one lower-triangular test matrix P. The objective is evaluated as
% kappa_inf(P*Q).
P = tril(randn(7,7));

% Regenerate P if it is numerically singular.
while rcond(P) < eps
    P = tril(randn(3,3));
end
K = size(P,1);

% Store figures and the final solutions under a dimension-specific folder.
plotDir = fullfile('plots',  sprintf('K=%d', K));

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% Each cell stores the output struct returned by givens_algorithm.
solutions = cell(N,1);

% Run independent searches in parallel. Each iteration starts from a new
% Latin-hypercube design inside givens_algorithm.
parfor iP = 1:N


    %% Problem setup

    nTheta = K*(K-1)/2;
    results = givens_algorithm(P, 10000, 100, @norm_infinity);
    solutions{iP} = results;
    %% Diagnostics

    kappa_identity = condition_number(P, @norm_infinity);
    orthogonality_error = norm(results.Q_best.'*results.Q_best - eye(K), 'fro');

    fprintf('\nGivens search completed for P %d.\n', iP)
    fprintf('K: %d\n', K)
    fprintf('Number of Givens angles: %d\n', nTheta)
    fprintf('kappa_inf(P): %.12g\n', kappa_identity)
    fprintf('best kappa_inf(P*Q): %.12g\n', results.kappa_best)
    fprintf('orthogonality error ||Q''Q-I||_F: %.3e\n', orthogonality_error)

    disp('theta_best =')
    disp(results.theta_best*(180/pi))

    disp('Q_best =')
    disp(results.Q_best)

    % For K=2, the objective can be shown as a one-dimensional curve. For
    % K>2, only the first three Givens angles are visualized.
    if nTheta==1
        figure
        plot(results.thetaCandidates, results.kappaCandidates, 'LineWidth', 1.2)
        hold on
        plot(results.theta_best(1), results.kappa_best, 'ro', 'MarkerFaceColor', 'r')
        xlabel('\theta')
        ylabel('\kappa_\infty(PQ)')
        title(sprintf('Givens angle exploration, P %03d', iP))
        grid on
    else
        figure
        scatter3( ...
            results.thetaCandidates(:, 1), ...
            results.thetaCandidates(:, 2), ...
            results.thetaCandidates(:, 3), ...
            18, log10(results.kappaCandidates), 'filled')
        hold on
        scatter3(results.theta_best(1), results.theta_best(2), results.theta_best(3), ...
            90, log10(results.kappa_best), 'r', 'filled')
        xlabel('\theta_{12}')
        ylabel('\theta_{13}')
        zlabel('\theta_{23}')
        title(sprintf('Givens angle exploration, P %03d', iP))
        cb = colorbar;
        ylabel(cb, 'log_{10}(\kappa_\infty(PQ))')
        grid on
        view(45, 25)
    end

    % The plotting code is kept for inspection, but saving each figure can
    % be expensive for large N.
    %saveas(gcf, fullfile(plotDir, sprintf('givens_K%d_P_%03d.png', K, iP)));
    close(gcf)

end

elapsed_time = toc;
fprintf('\nTotal elapsed time: %.3f seconds\n', elapsed_time)

% Save the complete set of solutions and the matrix P used in the search.
save(fullfile(plotDir,  sprintf('givens_K%d_solutions.mat', K)), 'solutions', 'P')
