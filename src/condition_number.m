function kappa_infinity = condition_number(A)

kappa_infinity= norm_infinity(A)*norm_infinity(inv(A));

end