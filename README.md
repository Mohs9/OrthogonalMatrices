# Orthogonal Matrices Optimization

This project develops algorithmic ideas for optimization problems constrained to
the space of orthogonal matrices.

The main objective is to study projected gradient and subgradient methods for
finding an orthogonal matrix \(Q\) that improves a structural decomposition while
preserving the reduced-form covariance matrix. In particular, the project
focuses on minimizing the infinity-norm condition number of \(PQ\), where \(P\)
is a Cholesky factor and \(Q\) belongs to the orthogonal group.

The optimization problem is written as

$$
Q^* = \text{argmin}_{Q \in \mathcal{O}(K)} \kappa_\infty(PQ).
$$

The proposed algorithm computes gradients or subgradients of the objective,
projects them onto the tangent space of the orthogonal group, and updates \(Q\)
along a Riemannian descent direction.

Project notes are available in `docs/algorithm_idea.qmd`, and MATLAB routines
are organized in `src/` and `scripts/`.
