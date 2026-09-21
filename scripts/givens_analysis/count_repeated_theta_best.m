% Count repeated theta_best points and Q_best matrices across solutions.
%
% The Givens parametrization is not unique: two different theta_best vectors
% may produce the same, or nearly the same, orthogonal matrix Q_best. This
% script reports repetitions both in angle space and in matrix space.

activate

% Load the stored solutions for the selected dimension.
K = 3;
plotDir = fullfile('plots', sprintf('K=%d', K));
load(fullfile(plotDir, sprintf('givens_K%d_solutions.mat', K)), 'solutions')

% Build one row per solution for both theta_best and vectorized Q_best.
N = numel(solutions);
nTheta = numel(solutions{1}.theta_best);
Theta_best = zeros(N, nTheta);
Q_best_vec = zeros(N, K*K);

for i = 1:N
    Theta_best(i,:) = wrap_to_pi_local(solutions{i}.theta_best(:).');
    Q_best_vec(i,:) = reshape(solutions{i}.Q_best, 1, []);
end

% Main tolerance used for the detailed repetition tables.
tol = 1e-4;
Theta_round = round(Theta_best/tol)*tol;
Q_round = round(Q_best_vec/tol)*tol;

% Count unique angle vectors after wrapping and rounding.
[Theta_unique, ~, idx_theta_group] = unique(Theta_round, 'rows');
theta_counts = accumarray(idx_theta_group, 1);

% Count unique matrices after vectorizing and rounding.
[Q_unique, ~, idx_Q_group] = unique(Q_round, 'rows');
Q_counts = accumarray(idx_Q_group, 1);

n_theta_unique = size(Theta_unique, 1);
n_theta_repeated = N - n_theta_unique;
n_Q_unique = size(Q_unique, 1);
n_Q_repeated = N - n_Q_unique;

% Sensitivity check: repeat the Q-counting exercise across several
% tolerances to see whether repetitions are exact or only approximate.
tols = [1e-6 1e-5 1e-4 1e-3 1e-2];

for t = tols
    Q_round = round(Q_best_vec/t)*t;
    n_Q_unique = size(unique(Q_round,'rows'),1);

    fprintf('tol = %.0e | unique Q = %d | repeated = %d\n', ...
        t, n_Q_unique, N - n_Q_unique)
end

fprintf('\nRepeated theta_best count\n')
fprintf('K: %d\n', K)
fprintf('Tolerance: %.1e\n', tol)
fprintf('Total solutions: %d\n', N)
fprintf('Unique theta_best points: %d\n', n_theta_unique)
fprintf('Repeated theta_best points: %d\n', n_theta_repeated)
fprintf('Unique Q_best matrices: %d\n', n_Q_unique)
fprintf('Repeated Q_best matrices: %d\n\n', n_Q_repeated)

% Tables store the group id and the number of solutions in each group.
theta_group = (1:n_theta_unique).';
theta_repeat_table = table(theta_group, theta_counts);

Q_group = (1:n_Q_unique).';
Q_repeat_table = table(Q_group, Q_counts);

disp('Theta groups with repetitions:')
disp(theta_repeat_table(theta_counts > 1,:))

disp('Q groups with repetitions:')
disp(Q_repeat_table(Q_counts > 1,:))

% Save MATLAB data for later inspection.
save(fullfile(plotDir, 'theta_best_repetition_count.mat'), ...
    'Theta_best', ...
    'Theta_round', ...
    'Theta_unique', ...
    'idx_theta_group', ...
    'theta_counts', ...
    'theta_repeat_table', ...
    'Q_best_vec', ...
    'Q_round', ...
    'Q_unique', ...
    'idx_Q_group', ...
    'Q_counts', ...
    'Q_repeat_table', ...
    'tol')

% Save compact CSV tables for quick review.
writetable(theta_repeat_table, ...
    fullfile(plotDir, 'theta_best_repetition_count.csv'))

writetable(Q_repeat_table, ...
    fullfile(plotDir, 'Q_best_repetition_count.csv'))
