function g = norm_infinity(A)

rowSumA = sum(abs(A),2);
g = max(rowSumA);


end