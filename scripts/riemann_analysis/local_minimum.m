
% Run a local Riemannian descent experiment from one random initial rotation.

% Activate the project paths.
activate

% Set the matrix dimension for this small example.
K = 2;

% Define a test matrix for the condition-number objective.
L = [5, 0;1,6];

% Generate a random 2-by-2 orthogonal rotation matrix.
theta = 2*pi*rand;

Q_0 = [cos(theta), -sin(theta);
    sin(theta),  cos(theta)];

% Set optimization parameters.
maxIter = 200;
tol = 1e-8;
alpha = 0.1;

% Compare the initial condition number with the optimized result.
kappa_0 = condition_number(L*Q_0, @norm_infinity);
[Q_ast, kappa] = riemman_algorithm(L, Q_0, alpha, maxIter, tol, @subgrad_norm_inf,@norm_infinity);
