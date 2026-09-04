
% 1. Definir el tamaño de la matriz cuadrada (por ejemplo, 4x4)
K = 2;

% 2. Generar una matriz aleatoria Q inicial completa

% Matriz ortogonal
L = [5, 0;1,4];
theta = 2*pi*rand;

Q_0 = [cos(theta), -sin(theta);
    sin(theta),  cos(theta)];

maxIter = 200;
tol = 1e-8;
alpha = 0.1;

kappa_0 = condition_number(L*Q_0);
[Q_ast, kappa] = riemman_algorithm(L, Q_0, alpha, maxIter, tol);

