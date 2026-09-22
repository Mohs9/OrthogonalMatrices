function [Q0, theta_best, kappa_best] = multistart_Q0(L, nStart)

K = size(L,1);
nTheta = K*(K-1)/2;

lb = -pi*ones(1,nTheta);
ub =  pi*ones(1,nTheta);

objective = @(theta) cond(L*construct_Givens_matrices(theta,K), inf);

options = optimoptions("fmincon", ...
    "Display","off", ...
    "Algorithm","sqp");

kappa_best = Inf;
theta_best = [];

for s = 1:nStart

    % Random initial angles
    theta_init = lb + (ub-lb).*rand(1,nTheta);

    [theta_sol,kappa_sol] = fmincon( ...
        objective, theta_init, ...
        [],[],[],[], ...
        lb,ub,[],options);

    if kappa_sol < kappa_best
        kappa_best = kappa_sol;
        theta_best = theta_sol;
    end
end

Q0 = construct_Givens_matrices(theta_best,K);

end


