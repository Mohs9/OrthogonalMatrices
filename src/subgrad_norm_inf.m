function [Gphi, phi, kappa] = subgrad_norm_inf(L,Q)

K = size(L,1);
Linv = L \ eye(K);

%% ||LQ||_1
A = L*Q;


colSumA = sum(abs(A),1);
g = max(colSumA);

tolA = 1e-12*max(1,g);
activeA = find(abs(colSumA-g) <= tolA);

EA = zeros(K);

for j = activeA.'
    EA(:,j) = sign(A(:,j))/length(activeA);
end

Gg = L.'*EA;


%% ||Q' L^{-1}||_1
B = Q.'*Linv;


colSumB = sum(abs(B),1);
h = max(colSumB);

tolB = 1e-12*max(1,h);
activeB = find(abs(colSumB-h) <= tolB);

EB = zeros(K);

for j = activeB.'
    EB(:,j) = sign(B(:,j))/length(activeB);
end

Gh = Linv*EB.';

%% Objective
kappa = g*h;
phi = log(g) + log(h);

%% Euclidean subgradient of log(kappa)
Gphi = Gg/g + Gh/h;

end