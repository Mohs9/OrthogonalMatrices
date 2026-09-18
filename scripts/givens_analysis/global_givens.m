%---------------------------------------------------------------------
% This script whe are going to explore all space of ortogonal matrices
% using givens rotation matrices.
%---------------------------------------------------------------------

activate
%rng(1, 'twister')
tic

N = 400;
P = tril(randn(4,4));

while rcond(P) < eps
    P = tril(randn(3,3));
end
K = size(P,1);

plotDir = fullfile('plots',  sprintf('K=%d', K));

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end


solutions = cell(N,1);
parfor iP = 1:N


    %% Problem setup


    nTheta = K*(K-1)/2;
    results = givens_algorithm(P, 10000, 100);
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
    %saveas(gcf, fullfile(plotDir, sprintf('givens_K%d_P_%03d.png', K, iP)));
    close(gcf)

end

elapsed_time = toc;
fprintf('\nTotal elapsed time: %.3f seconds\n', elapsed_time)

save(fullfile(plotDir,  sprintf('givens_K%d_solutions.mat', K)), 'solutions', 'P')

