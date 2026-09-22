function [value, grad, kappa, Omega] = objective_smooth(Q, P, norm_func)
%OBJECTIVE_SMOOTH Log objective and Riemannian subgradient.
%
%   Reuse the exact subgradient already implemented for the Riemannian
%   algorithm. The name is kept so riemannian_refine stays simple.

if nargin < 3 || isempty(norm_func)
    norm_func = @norm_1;
end

normName = func2str(norm_func);

if strcmp(normName, 'norm_1')
    [Gphi, value, kappa] = subgrad_norm_inf(P, Q);
elseif strcmp(normName, 'norm_infinity')
    [Gphi, value, kappa] = subgrad_norm_1(P, Q);
else
    error('Riemannian refinement only supports @norm_1 and @norm_infinity.');
end

Omega = skew(Q' * Gphi);
grad = Q * Omega;
end
