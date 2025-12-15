function setup_project()
%SETUP_PROJECT Adds project directories to the MATLAB path.
% 
%   Run this script once before executing the master harness or any models.

    % Get root directory of the current script
    root_dir = fileparts(mfilename('fullpath'));
    
    % Define paths to add
    paths_to_add = {
        fullfile(root_dir, 'api'),
        fullfile(root_dir, 'models', 'itm'),
        fullfile(root_dir, 'models', 'empirical'),
        fullfile(root_dir, 'models', 'diffraction'),
        fullfile(root_dir, 'models', 'fullwave'),
        fullfile(root_dir, 'data'),
        fullfile(root_dir, 'output')
    };
    
    % Add paths
    fprintf('Adding project paths...\n');
    for i = 1:length(paths_to_add)
        p = paths_to_add{i};
        if exist(p, 'dir')
            addpath(p);
            fprintf('  + %s\n', p);
        else
            fprintf('  ! Warning: Directory not found: %s\n', p);
        end
    end
    
    savepath;
    fprintf('Path setup complete.\n');

end
