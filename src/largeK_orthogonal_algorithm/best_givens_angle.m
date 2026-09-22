function [theta_best, f_best] = best_givens_angle(Q, P, R, i, j, norm_func, opts)
K = size(Q, 1);
nGrid = opts.GIVENS_GRID_SIZE;

if mod(nGrid, 2) == 0
    nGrid = nGrid + 1;
end

theta_grid = linspace(-opts.GIVENS_RADIUS, opts.GIVENS_RADIUS, nGrid);
f_values = zeros(nGrid, 1);

for k = 1:nGrid
    Qtrial = Q * local_givens_rotation(K, i, j, theta_grid(k));
    f_values(k) = objective_exact(Qtrial, P, R, norm_func);
end

[f_best, idx_best] = min(f_values);
theta_best = theta_grid(idx_best);
end
