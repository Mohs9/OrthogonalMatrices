function [Q_new, kappa_new]= riemman_algorithm(L,Q_0,alpha, maxIter, tol, subgrad_norm, norm_p)
%RIEMMAN_ALGORITHM Run a Riemannian descent method on the orthogonal group.
%   Starting from Q_0, the algorithm computes a Euclidean subgradient,
%   projects it onto the tangent space, and updates Q by a matrix exponential.
%   The step size is adapted by accepting only condition-number improvements.

% Initialize the current iterate and line-search bounds.
Q_t = Q_0;
alpha_max = 0.1;
alpha_min = 1e-10;

for iter = 1:maxIter

    % Evaluate the Euclidean subgradient and objective at the current point.
    [Gphi, ~, kappa_t] = subgrad_norm(L,Q_t);

    % Build the skew-symmetric tangent generator Omega.
    Omega = skew( Gphi*Q_t');

    % Stop when the tangent generator is sufficiently small.
    if norm(Omega,'fro') < tol
        break
    end

    % Move along the manifold using the exponential map.
    Q_new =  expm(-alpha * Omega) * Q_t ;

    % Evaluate the condition number after the candidate update.
    kappa_new = condition_number(L * Q_new, norm_p);

    if kappa_new < kappa_t

        % Accept the update and increase the step size slightly.
        Q_t = Q_new;
        alpha = min(1.2*alpha, alpha_max);

    else

        % Reject the update and reduce the step size.
        alpha = alpha/2;

        % Stop if the step size has become too small to be useful.
        if alpha < alpha_min
            break
        end

    end

end
end
