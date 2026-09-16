function kappa_p = condition_number(A, norm_p)
%CONDITION_NUMBER Compute the p-norm condition number of a matrix.
%   The value is ||A||_p * ||A^{-1}||_p, using the project-specific
%   p norm implementation.

kappa_p = norm_p(A)*norm_p(inv(A));

end
