% ---------------------------------------------------------------------
% Global multistart experiment for the Riemannian descent algorithm.
%
% The script draws many random orthogonal initial matrices Q_0, runs the
% Riemannian algorithm from each starting point, and keeps the best local
% solution found for the infinity-norm condition number objective.
% ---------------------------------------------------------------------

activate;
tic;

N = 10000; % Number of orthogonal initial matrices to evaluate.
m = 9;    % Number of rows.
n = 9;    % Number of columns, with m >= n.

% Store the optimized matrices and their objective values.
Q_resultados = cell(N, 1);
f_resultados = zeros(N, 1);

% Set optimization parameters.
maxIter = 10000;
tol = 1e-8;
alpha = 0.2;

% Generate a random test matrix for the condition-number objective.
L = randn(9,9);

% Run the local Riemannian descent from many random initial conditions.
parfor i = 1:N
    % Generate a random orthogonal initial matrix.
    Q_0 = orthogonal_matrix_generator(m, n);

    % Minimize the objective from this initial matrix.
    [Q_ast, kappa] = riemman_algorithm(L, Q_0, alpha, maxIter, tol, @subgrad_norm_inf, @norm_infinity);

    % Store the local solution and objective value.
    Q_resultados{i} = Q_ast;
    f_resultados(i) = kappa;
end

% Select the best local result among all restarts.
[f_global, idx_mejor] = min(f_resultados);
Q_global = Q_resultados{idx_mejor};

time = toc;
