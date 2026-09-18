% Plot the objective over the first three Givens angles.

activate

K = 5;
plotDir = fullfile('plots', sprintf('K=%d', K));
load(fullfile(plotDir, sprintf('givens_K%d_solutions.mat', K)), 'solutions')

Theta = [];
kappa = [];

N = numel(solutions);
Theta_best = zeros(N, K);
kappa_best = zeros(N, 1);

for i = 1:N
    Theta = [Theta; solutions{i}.thetaCandidates(:,1:K)];
    kappa = [kappa; solutions{i}.kappaCandidates];

    Theta_best(i,:) = solutions{i}.theta_best(1:K);
    kappa_best(i) = solutions{i}.kappa_best;
end

figure
%scatter3(Theta(:,1), Theta(:,2), Theta(:,3), ...
% 8, log10(kappa), 'filled')
hold on
scatter3(Theta_best(:,3), Theta_best(:,4), Theta_best(:,5), ...
    80, 'r', 'filled')

xlabel('\theta_1')
ylabel('\theta_2')
zlabel('\theta_3')
title(sprintf('Objective over first three Givens angles, K=%d', K))
cb = colorbar;
ylabel(cb, 'log_{10}(\kappa_\infty(PQ))')
grid on
view(45, 25)

saveas(gcf, fullfile(plotDir, 'objective_last3_angles.png'))
