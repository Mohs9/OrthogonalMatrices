function matrix = skew(A)
%SKEW Return the skew-symmetric part of a square matrix.

matrix = (A-A')/2;
end
