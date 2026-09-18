% ---------------------------------------------------------------------
% Check whether Q_best solutions are the same up to signed permutations.
% ---------------------------------------------------------------------

clearvars -except solutions
activate;
plotDir = fullfile('plots', 'K=3');
solutionsFile = fullfile(plotDir, 'givens_K3_solutions.mat');

if ~exist('solutions', 'var')
    load(solutionsFile, 'solutions')
end

N = numel(solutions);
distance_perm = zeros(N);

for i = 1:N
    Q_i = solutions{i}.Q_best;

    for j = i+1:N
        Q_j = solutions{j}.Q_best;

        distance_perm(i,j) = nearest_signed_permutation(Q_i, Q_j);
        distance_perm(j,i) = distance_perm(i,j);
    end
end

tol = 1;
same_perm = distance_perm <= tol;

nPairs = N*(N-1)/2;
nSamePairs = nnz(triu(same_perm, 1));

fprintf('\nSigned permutation check\n')
fprintf('Equivalent pairs: %d of %d\n', nSamePairs, nPairs)
fprintf('Share equivalent: %.2f%%\n', 100*nSamePairs/nPairs)
fprintf('Max distance: %.3e\n', max(distance_perm(:)))

figure
imagesc(distance_perm)
axis image
colorbar
xlabel('solution id')
ylabel('solution id')
title('Distance after best signed permutation')
saveas(gcf, fullfile(plotDir, 'signed_permutation_distance_matrix.png'))

figure
imagesc(same_perm)
axis image
colormap(gray)
xlabel('solution id')
ylabel('solution id')
title(sprintf('Same up to signed permutation, tol %.1e', tol))
saveas(gcf, fullfile(plotDir, 'signed_permutation_equivalence_matrix.png'))


pairDistances = distance_perm(triu(true(N), 1));
kappas = cellfun(@(s) s.kappa_best, solutions);

figure
histogram(pairDistances, 30)
xlabel('min_S ||Q_j - Q_i S||_F')
ylabel('count')
title('Pairwise distances after best signed permutation')
grid on
saveas(gcf, fullfile(plotDir, 'signed_permutation_distance_histogram.png'))

figure
histogram(kappas, 30)
xlabel('\kappa_\infty(PQ_{best})')
ylabel('count')
title('Distribution of best condition numbers')
grid on
saveas(gcf, fullfile(plotDir, 'kappa_best_histogram.png'))


%%%%------------------------- compare Givens angle -------

nTheta = numel(solutions{1}.theta_best);
Theta = zeros(N, nTheta);

for i = 1:N
    Theta(i,:) = wrap_to_pi_local(solutions{i}.theta_best(:).');
end

theta_tol = 1e-6;
common_angle_result = detect_common_angle(Theta, theta_tol);

fprintf('\nCommon angle detection\n')
fprintf('alpha: %.12g\n', common_angle_result.alpha)
fprintf('alpha/pi: %.12g\n', common_angle_result.alpha_over_pi)
fprintf('RMS residual: %.3e\n', common_angle_result.rmsResidual)
fprintf('Max residual: %.3e\n', common_angle_result.maxResidual)
fprintf('Passes tolerance: %d\n', common_angle_result.passesTolerance)
fprintf('Selection rule: %s\n', common_angle_result.selectionRule)

save(fullfile(plotDir, 'common_angle_detection.mat'), ...
    'Theta', ...
    'theta_tol', ...
    'common_angle_result')

givens_angle_results = compare_Givens_angles( ...
    Theta, theta_tol, common_angle_result.alpha);

fprintf('\nGivens angle comparison\n')
fprintf('Number of Givens angles: %d\n', nTheta)
fprintf('Equivalent angle pairs modulo detected alpha: %d of %d\n', ...
    height(givens_angle_results.pairs), nPairs)

save(fullfile(plotDir, 'givens_angle_comparison.mat'), ...
    'Theta', ...
    'theta_tol', ...
    'givens_angle_results')

figure
imagesc(givens_angle_results.distance)
axis image
colorbar
xlabel('solution id')
ylabel('solution id')
title('Givens angle distance')
saveas(gcf, fullfile(plotDir, 'givens_angle_distance_matrix.png'))

figure
imagesc(givens_angle_results.equivalentPi2)
axis image
colormap(gray)
xlabel('solution id')
ylabel('solution id')
title(sprintf('Givens angles equivalent modulo %.4g*pi, tol %.1e', ...
    common_angle_result.alpha_over_pi, theta_tol))
saveas(gcf, fullfile(plotDir, 'givens_angle_equivalence_pi2_matrix.png'))


figure
plot(common_angle_result.alphaGrid/pi, ...
    common_angle_result.maxResidualGrid, 'LineWidth', 1.2)
hold on
plot(common_angle_result.alpha_over_pi, ...
    common_angle_result.maxResidual, 'ro', 'MarkerFaceColor', 'r')
yline(theta_tol, '--k', 'LineWidth', 1.0)
xlabel('\alpha/\pi')
ylabel('maximum residual')
title('Common angle detection')
grid on
saveas(gcf, fullfile(plotDir, 'common_angle_detection_residuals.png'))
