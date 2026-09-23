% Analyze the optimized kappa_star values across stored Haar + Givens solutions.

activate

% Load the stored solutions for the selected dimension.
K = 25;
plotDir = fullfile('plots', sprintf('K=%d', K));
load(fullfile(plotDir, sprintf('largeK_orthogonal_K%d_solutions.mat', K)))

N = numel(solutions);
kappas = cellfun(@(s) s.kappa_star, solutions);
eigenvalues_Q = complex(zeros(N, K));

% Tolerance for grouping numerically equal objective values.
tol = 1e-4;
kappas_round = round(kappas/tol)*tol;

% Print Q_star matrices for groups with repeated kappa_star values. Set this
% to true to print singleton groups as well.
display_singletons = false;

for i = 1:N
    eigenvalues_Q(i,:) = eig(solutions{i}.Q_star).';
end

eigenvalue_moduli = abs(eigenvalues_Q);
eigenvalue_angles = angle(eigenvalues_Q);

[unique_kappas, ~, idx_kappa_group] = unique(kappas_round);
kappa_counts = accumarray(idx_kappa_group, 1);
nGroups = numel(unique_kappas);
group_solution_indices = cell(nGroups, 1);
Q_star_groups = cell(nGroups, 1);

for iGroup = 1:nGroups
    group_solution_indices{iGroup} = find(idx_kappa_group == iGroup);
    group_indices = group_solution_indices{iGroup};
    Q_group = zeros(K, K, numel(group_indices));

    for iMember = 1:numel(group_indices)
        Q_group(:,:,iMember) = solutions{group_indices(iMember)}.Q_star;
    end

    Q_star_groups{iGroup} = Q_group;
end

fprintf('\nkappa_star analysis\n')
fprintf('K: %d\n', K)
fprintf('Total solutions: %d\n', N)
fprintf('Tolerance: %.1e\n', tol)
fprintf('min(kappa_star): %.15g\n', min(kappas))
fprintf('max(kappa_star): %.15g\n', max(kappas))
fprintf('range: %.3e\n', max(kappas)-min(kappas))
fprintf('mean(kappa_star): %.15g\n', mean(kappas))
fprintf('std(kappa_star): %.3e\n', std(kappas))
fprintf('Unique kappa groups: %d\n', nGroups)
fprintf('max |abs(eig(Q_star))-1|: %.3e\n', ...
    max(abs(eigenvalue_moduli(:)-1)))

kappa_table = table((1:nGroups).', unique_kappas, kappa_counts, ...
    'VariableNames', {'Group', 'KappaStar', 'Count'});

disp(kappa_table)

fprintf('\nQ_star matrices by kappa_star group\n')

printedGroups = 0;
group_id_detail = [];
solution_id_detail = [];
kappa_detail = [];

for iGroup = 1:nGroups
    group_indices = group_solution_indices{iGroup};

    group_id_detail = [group_id_detail; repmat(iGroup, numel(group_indices), 1)]; %#ok<AGROW>
    solution_id_detail = [solution_id_detail; group_indices(:)]; %#ok<AGROW>
    kappa_detail = [kappa_detail; kappas(group_indices(:))]; %#ok<AGROW>

    if ~display_singletons && numel(group_indices) == 1
        continue
    end

    printedGroups = printedGroups + 1;
    fprintf('\nGroup %d\n', iGroup)
    fprintf('Rounded kappa_star: %.15g\n', unique_kappas(iGroup))
    fprintf('Solutions in group: ')
    fprintf('%d ', group_indices)
    fprintf('\n')

    for iMember = 1:numel(group_indices)
        iSolution = group_indices(iMember);

        fprintf('\nSolution %d, kappa_star = %.15g\n', ...
            iSolution, kappas(iSolution))
        disp('Q_star =')
        disp(solutions{iSolution}.Q_star)
    end
end

if printedGroups == 0
    fprintf('No groups with more than one solution were found. ')
    fprintf('Set display_singletons = true to print all Q_star matrices.\n')
end

group_detail_table = table(group_id_detail, solution_id_detail, kappa_detail, ...
    'VariableNames', {'Group', 'Solution', 'KappaStar'});

solution_id = repelem((1:N).', K);
eigenvalue_id = repmat((1:K).', N, 1);
eigenvalue_real = real(eigenvalues_Q(:));
eigenvalue_imag = imag(eigenvalues_Q(:));
eigenvalue_modulus = eigenvalue_moduli(:);
eigenvalue_angle = eigenvalue_angles(:);

eigenvalue_table = table( ...
    solution_id, ...
    eigenvalue_id, ...
    eigenvalue_real, ...
    eigenvalue_imag, ...
    eigenvalue_modulus, ...
    eigenvalue_angle, ...
    'VariableNames', { ...
    'Solution', ...
    'Eigenvalue', ...
    'Real', ...
    'Imag', ...
    'Modulus', ...
    'Angle'});

save(fullfile(plotDir, 'kappa_star_analysis.mat'), ...
    'K', ...
    'kappas', ...
    'tol', ...
    'kappas_round', ...
    'unique_kappas', ...
    'idx_kappa_group', ...
    'kappa_counts', ...
    'kappa_table', ...
    'group_solution_indices', ...
    'Q_star_groups', ...
    'group_detail_table', ...
    'display_singletons', ...
    'eigenvalues_Q', ...
    'eigenvalue_moduli', ...
    'eigenvalue_angles', ...
    'eigenvalue_table')

writetable(kappa_table, fullfile(plotDir, 'kappa_star_analysis.csv'))
writetable(group_detail_table, fullfile(plotDir, 'kappa_star_group_members.csv'))
writetable(eigenvalue_table, fullfile(plotDir, 'Q_star_eigenvalues.csv'))
