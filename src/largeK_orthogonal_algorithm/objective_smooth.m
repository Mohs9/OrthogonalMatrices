function [value, grad, kappa, Omega] = objective_smooth(Q, P, norm_func)
%OBJECTIVE_SMOOTH Log objective and Riemannian subgradient.
%
%   Returns value = log(kappa), a projected subgradient grad on the tangent
%   space of O(K), the condition number kappa, and the skew-symmetric matrix
%   Omega used to update Q with a matrix exponential.

if nargin < 3 || isempty(norm_func)
    norm_func = @norm_1;
end

normName = func2str(norm_func);

% The helper routines are named after the dual-term subgradient, so @norm_1
% uses subgrad_norm_inf and vice versa.
if strcmp(normName, 'norm_1')
    [Gphi, value, kappa] = subgrad_norm_inf(P, Q);
elseif strcmp(normName, 'norm_infinity')
    [Gphi, value, kappa] = subgrad_norm_1(P, Q);
else
    error('Riemannian refinement only supports @norm_1 and @norm_infinity.');
end

% Project the Euclidean subgradient Gphi onto the tangent space at Q.
Omega = skew(Q' * Gphi);
grad = Q * Omega;
end
