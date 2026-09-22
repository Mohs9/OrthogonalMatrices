function [kappa, ell] = objective_exact(Q, P, R, norm_func)
%OBJECTIVE_EXACT Exact kappa_1(P*Q), avoiding inv(P*Q).

norm_A = norm_func(P*Q);
norm_C = norm_func(Q'*R);
kappa = norm_A * norm_C;
ell = log(norm_A) + log(norm_C);
end
