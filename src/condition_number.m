function kappa_infinity = condition_number(A)
%CONDITION_NUMBER Compute the infinity-norm condition number of a matrix.
%   The value is ||A||_inf * ||A^{-1}||_inf, using the project-specific
%   infinity norm implementation.

kappa_infinity = norm_infinity(A)*norm_infinity(inv(A));

end
