% ---------------------------------------------------------------------
% Setup for a K=4 Givens search.
%
% This script currently prepares the matrix P, output directory, and storage
% cell array for a K=4 experiment. The optimization loop is not included in
% the current file.
% ---------------------------------------------------------------------

activate
%rng(1, 'twister')
tic

% Number of independent searches intended for this experiment.
N = 400;

% Generate one lower-triangular matrix P for the condition-number objective.
P = tril(randn(4,4));

% Regenerate P if it is numerically singular.
while rcond(P) < eps
    P = tril(randn(3,3));
end
K = size(P,1);

% Store outputs under a folder named after the matrix dimension.
plotDir = fullfile('plots',  sprintf('K=%d', K));

if ~exist(plotDir, 'dir')
    mkdir(plotDir);
end

% Placeholder for the output structs returned by givens_algorithm.
solutions = cell(N,1);


norm_4(P, 1000)