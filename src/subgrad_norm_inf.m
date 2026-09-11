function [Gphi, phi, kappa] = subgrad_norm_inf(L,Q)
%SUBGRAD_NORM_INF Compute a Euclidean subgradient of log kappa for column sums.
%   The routine evaluates g = ||LQ||_1 and h = ||Q'L^{-1}||_1, then
%   returns a subgradient of phi = log(g) + log(h) with respect to Q.

K = size(L,1);
Linv = L \ eye(K);

%% ||LQ||_1
A = L*Q;

% Find the active columns that attain the 1-norm.
colSumA = sum(abs(A),1);
g = max(colSumA);

tolA = 1e-12*max(1,g);
activeA = find(abs(colSumA-g) <= tolA);

% Average the column-wise sign patterns across active columns.
EA = zeros(K);

for j = activeA.'
    EA(:,j) = sign(A(:,j))/length(activeA);
end

% Chain rule for A = LQ.
Gg = L.'*EA;


%% ||Q' L^{-1}||_1
B = Q.'*Linv;

% Find the active columns of the inverse term.
colSumB = sum(abs(B),1);
h = max(colSumB);

tolB = 1e-12*max(1,h);
activeB = find(abs(colSumB-h) <= tolB);

% Average the column-wise sign patterns across active columns.
EB = zeros(K);

for j = activeB.'
    EB(:,j) = sign(B(:,j))/length(activeB);
end

% Chain rule for B = Q'L^{-1}.
Gh = Linv*EB.';

%% Objective
kappa = g*h;
phi = log(g) + log(h);

%% Euclidean subgradient of log(kappa)
Gphi = Gg/g + Gh/h;

end
