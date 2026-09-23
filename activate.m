% Create an environment variable to hold the project root path
OrthogonalMatrices_PROJECTDIR = pwd();
setenv("OrthogonalMatrices_PROJECTDIR", OrthogonalMatrices_PROJECTDIR);
disp(['Project root folder: ', OrthogonalMatrices_PROJECTDIR]);

% Add src folder with subfolders to the path
addpath(genpath(fullfile(OrthogonalMatrices_PROJECTDIR, 'src')));

