% Create an environment variable to hold the project root path
OrthogonalMatrices_PROJECTDIR = pwd();
setenv("OrthogonalMatrices_PROJECTDIR", OrthogonalMatrices_PROJECTDIR);
disp(['Project root folder: ', OrthogonalMatrices_PROJECTDIR]);

% Add src folder with subfolders to the path
addpath(genpath(fullfile(OrthogonalMatrices_PROJECTDIR, 'src')));

% Trying to initialize the IRIS toolbox.
try
    irisstartup;
catch ME
    disp('Perhaps you need to add the IRIS toolbox to the path...');
    rethrow(ME);
end

% Current version of the IRIS toolbox used in the project.
IRIS_VERSION = '20191112';

if ~strcmp(IRIS_VERSION, irisversion)
    warning(...
        sprintf('IRIS toolbox version %s is recommended for this project.', IRIS_VERSION) ...
        )
end

% Configuring the default time series constructor to tseries
iris.set('DefaultTimeSeriesConstructor', @tseries)