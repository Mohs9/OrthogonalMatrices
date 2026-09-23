function [theta_best, f_best] = best_givens_angle(Q, P, R, i, j, norm_func, opts)
%BEST_GIVENS_ANGLE Find the best local angle in one Givens plane.
%
%   Evaluates the exact objective on a symmetric angle grid around zero and
%   returns the angle that gives the smallest value. Q is not modified here;
%   exact_givens_polish applies the winning rotation.

K = size(Q, 1);
nGrid = opts.GIVENS_GRID_SIZE;

% An odd number of grid points guarantees that theta = 0 is included.
if mod(nGrid, 2) == 0
    nGrid = nGrid + 1;
end

theta_grid = linspace(-opts.GIVENS_RADIUS, opts.GIVENS_RADIUS, nGrid);
f_values = zeros(nGrid, 1);

% Test rotations Q*G(i,j,theta) around the current solution.
for k = 1:nGrid
    Qtrial = Q * local_givens_rotation(K, i, j, theta_grid(k));
    f_values(k) = objective_exact(Qtrial, P, R, norm_func);
end

% The grid minimum defines the local candidate for this plane.
[f_best, idx_best] = min(f_values);
theta_best = theta_grid(idx_best);
end
