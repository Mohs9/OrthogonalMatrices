function opts = fill_default_options_largeK(opts)
%FILL_DEFAULT_OPTIONS_LARGEK Add defaults without overwriting user choices.

if ~isfield(opts, 'N_HAAR') || isempty(opts.N_HAAR)
    opts.N_HAAR = 5000;
end

if ~isfield(opts, 'N_ELITE') || isempty(opts.N_ELITE)
    opts.N_ELITE = 10;
end

if ~isfield(opts, 'SEED')
    opts.SEED = [];
end

if ~isfield(opts, 'RIEMANN_MAX_ITERS') || isempty(opts.RIEMANN_MAX_ITERS)
    opts.RIEMANN_MAX_ITERS = 100;
end

if ~isfield(opts, 'RIEMANN_GRAD_TOL') || isempty(opts.RIEMANN_GRAD_TOL)
    opts.RIEMANN_GRAD_TOL = 1e-7;
end

if ~isfield(opts, 'RIEMANN_REL_TOL') || isempty(opts.RIEMANN_REL_TOL)
    opts.RIEMANN_REL_TOL = 1e-10;
end

if ~isfield(opts, 'ALPHA0') || isempty(opts.ALPHA0)
    opts.ALPHA0 = 1;
end

if ~isfield(opts, 'LINESEARCH_RHO') || isempty(opts.LINESEARCH_RHO)
    opts.LINESEARCH_RHO = 0.5;
end

if ~isfield(opts, 'MAX_LINESEARCH') || isempty(opts.MAX_LINESEARCH)
    opts.MAX_LINESEARCH = 20;
end

if ~isfield(opts, 'GIVENS_MAX_SWEEPS') || isempty(opts.GIVENS_MAX_SWEEPS)
    opts.GIVENS_MAX_SWEEPS = 2;
end

if ~isfield(opts, 'GIVENS_RADIUS') || isempty(opts.GIVENS_RADIUS)
    opts.GIVENS_RADIUS = 0.05;
end

if ~isfield(opts, 'GIVENS_GRID_SIZE') || isempty(opts.GIVENS_GRID_SIZE)
    opts.GIVENS_GRID_SIZE = 9;
end

if ~isfield(opts, 'TOL_GIVENS') || isempty(opts.TOL_GIVENS)
    opts.TOL_GIVENS = 1e-10;
end

if ~isfield(opts, 'TOL_IMPROVEMENT') || isempty(opts.TOL_IMPROVEMENT)
    opts.TOL_IMPROVEMENT = 1e-12;
end
end
