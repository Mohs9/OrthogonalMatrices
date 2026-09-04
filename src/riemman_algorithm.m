function [Q_new, kappa_t]= riemman_algorithm(L,Q_0,alpha, maxIter, tol)

Q_t = Q_0;
alpha_max = 0.1;
alpha_min = 1e-10;

for iter = 1:maxIter

    [Gphi, ~, kappa_t] = subgrad_cond_inf(L,Q_t);

    Omega = skew(Q_t' * Gphi);

    if norm(Omega,'fro') < tol
        break
    end

    Q_new = Q_t * expm(-alpha * Omega);

    kappa_new = condition_number(L * Q_new);

    if kappa_new < kappa_t

        % Accept
        Q_t = Q_new;

        % Increase step size slightly
        alpha = min(1.2*alpha, alpha_max);

    else

        % Reject
        alpha = alpha/2;

        if alpha < alpha_min
            break
        end

    end

end
end