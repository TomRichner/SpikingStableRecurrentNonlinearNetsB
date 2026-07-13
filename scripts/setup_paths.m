function setup_paths()
%SETUP_PATHS Add the repository src (and scripts) directories to the MATLAB path.
%
% Call this helper before running any script or test that relies on code
% under src/ (sibling of scripts/). Mirrors the FractionalResevoir bootstrap.

scriptDir   = fileparts(mfilename('fullpath'));
projectRoot = fileparts(scriptDir);
srcPath     = fullfile(projectRoot, 'src');

if ~isfolder(srcPath)
    error('setup_paths:MissingSrc', ...
        'Could not find src directory at %s', srcPath);
end

addpath(genpath(scriptDir));
addpath(genpath(srcPath));
end
