function [Gphi, phi, kappa] = subgrad_cond_inf(L,Q)

K = size(L,1);
Linv = L \ eye(K);

%% ||LQ||_inf
A = L*Q;


rowSumA = sum(abs(A),2);
g = max(rowSumA);

tolA = 1e-12*max(1,g);
activeA = find(abs(rowSumA-g) <= tolA);

EA = zeros(K);

for i = activeA.'
    EA(i,:) = sign(A(i,:))/length(activeA);
end

Gg = L.'*EA;

%% ||Q' L^{-1}||_inf
B = Q.'*Linv;


rowSumB = sum(abs(B),2);
h = max(rowSumB);

tolB = 1e-12*max(1,h);
activeB = find(abs(rowSumB-h) <= tolB);

EB = zeros(K);

for i = activeB.'
    EB(i,:) = sign(B(i,:))/length(activeB);
end

Gh = Linv*EB.';

%% Objective
kappa = g*h;
phi = log(g) + log(h);

%% Euclidean subgradient of log(kappa)
Gphi = Gg/g + Gh/h;

end