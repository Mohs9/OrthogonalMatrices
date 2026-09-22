function opts = fill_default_options(opts)
%FILL_DEFAULT_OPTIONS Add defaults without overwriting user choices.

if ~isfield(opts, 'N_HAAR') || isempty(opts.N_HAAR)
    opts.N_HAAR = 50000;
end

if ~isfield(opts, 'N_ELITE') || isempty(opts.N_ELITE)
    opts.N_ELITE = 30;
end

if ~isfield(opts, 'DELTA0') || isempty(opts.DELTA0)
    opts.DELTA0 = 0.20;
end

if ~isfield(opts, 'TOL_STEP') || isempty(opts.TOL_STEP)
    opts.TOL_STEP = 1e-8;
end

if ~isfield(opts, 'TOL_IMPROVEMENT') || isempty(opts.TOL_IMPROVEMENT)
    opts.TOL_IMPROVEMENT = 1e-12;
end

if ~isfield(opts, 'MAX_SWEEPS') || isempty(opts.MAX_SWEEPS)
    opts.MAX_SWEEPS = 2000;
end

if ~isfield(opts, 'N_TANGENT') || isempty(opts.N_TANGENT)
    opts.N_TANGENT = 100;
end

if ~isfield(opts, 'TANGENT_RADIUS') || isempty(opts.TANGENT_RADIUS)
    opts.TANGENT_RADIUS = 0.05;
end

if ~isfield(opts, 'TOL_TANGENT') || isempty(opts.TOL_TANGENT)
    opts.TOL_TANGENT = 1e-8;
end

if ~isfield(opts, 'MAX_TANGENT_ROUNDS') || isempty(opts.MAX_TANGENT_ROUNDS)
    opts.MAX_TANGENT_ROUNDS = 20;
end

if ~isfield(opts, 'SEED')
    opts.SEED = [];
end
end
