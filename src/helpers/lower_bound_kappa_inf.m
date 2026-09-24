function [LB, LB_simple] = lower_bound_kappa_inf(L)
% LOWER_BOUND_KAPPA_INF
% Calcula cotas inferiores para kappa_inf(LQ),
% donde Q es cualquier matriz ortogonal.
%
% Inputs:
%   L : matriz cuadrada invertible n x n
%
% Outputs:
%   LB        : cota basada en normas de filas y Frobenius
%   LB_simple : cota simple kappa_2(L)/n

[n,m] = size(L);

if n ~= m
    error('L debe ser una matriz cuadrada.');
end

% Máxima norma 2 entre las filas de L
max_row_norm = max(vecnorm(L, 2, 2));

% L^{-1}, calculada resolviendo L*X = I
% Es preferible a usar inv(L)
Linv = L \ eye(n);

% Norma de Frobenius de L^{-1}
Linv_F = norm(Linv, 'fro');

% Cota inferior más fuerte
LB = max_row_norm * Linv_F / sqrt(n);

% Cota inferior simple: kappa_2(L)/n
LB_simple = cond(L, 2) / n;

end