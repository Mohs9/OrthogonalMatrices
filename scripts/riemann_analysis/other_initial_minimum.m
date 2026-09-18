
N = 1000;
% Store the optimized matrices and their objective values.
Q_resultados = cell(N, 1);
f_resultados = zeros(N, 1);

% Set optimization parameters.
maxIter = 100;
tol = 1e-8;
alpha = 0.2;

% Generate a random test matrix for the condition-number objective.
L = randn(9,9);

tic;
parfor i = 1:N
    % Generate a random orthogonal initial matrix.
    [Q_0, theta0, kappa0] =  multistart_Q0(L,10);

    % Minimize the objective from this initial matrix.
    [Q_ast, kappa] = riemman_algorithm(L, Q_0, alpha, maxIter, tol, @subgrad_norm_inf, @norm_infinity);

    % Store the local solution and objective value.
    Q_resultados{i} = Q_ast;
    f_resultados(i) = kappa;
end
tiempo_total = toc;
% Select the best local result among all restarts.
[f_global, idx_mejor] = min(f_resultados);
Q_global = Q_resultados{idx_mejor};