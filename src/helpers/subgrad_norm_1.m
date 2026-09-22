function [Gphi, phi, kappa] = subgrad_norm_1(L,Q)
%SUBGRAD_NORM_1 Compute a Euclidean subgradient of log kappa for row sums.
%   The routine evaluates g = ||LQ||_1 and h = ||Q'L^{-1}||_1, then
%   returns a subgradient of phi = log(g) + log(h) with respect to Q.

K = size(L,1);
Linv = L \ eye(K);

%% ||LQ||_inf
A = L*Q;

% Find the active rows that attain the infinity norm.
rowSumA = sum(abs(A),2);
g = max(rowSumA);

tolA = 1e-12*max(1,g);
activeA = find(abs(rowSumA-g) <= tolA);

% Average the row-wise sign patterns across active rows.
EA = zeros(K);

for i = activeA.'
    EA(i,:) = sign(A(i,:))/length(activeA);
end

% Chain rule for A = LQ.
Gg = L.'*EA;

%% ||Q' L^{-1}||_1
B = Q.'*Linv;

% Find the active rows of the inverse term.
rowSumB = sum(abs(B),2);
h = max(rowSumB);

tolB = 1e-12*max(1,h);
activeB = find(abs(rowSumB-h) <= tolB);

% Average the row-wise sign patterns across active rows.
EB = zeros(K);

for i = activeB.'
    EB(i,:) = sign(B(i,:))/length(activeB);
end

% Chain rule for B = Q'L^{-1}.
Gh = Linv*EB.';

%% Objective
kappa = g*h;
phi = log(g) + log(h);

%% Euclidean subgradient of log(kappa)
Gphi = Gg/g + Gh/h;

end
