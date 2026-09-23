function [kappa, ell] = objective_exact(Q, P, R, norm_func)
%OBJECTIVE_EXACT Evaluate kappa(P*Q) exactly without forming inv(P*Q).
%
%   With P such that Sigma_e = P*P' and R = P\I, we have
%   inv(P*Q) = Q'*R because Q is orthogonal. Therefore the condition number
%   is computed as ||P*Q|| * ||Q'*R||, and ell stores its logarithm.

norm_A = norm_func(P*Q);
norm_C = norm_func(Q'*R);
kappa = norm_A * norm_C;
ell = log(norm_A) + log(norm_C);
end
