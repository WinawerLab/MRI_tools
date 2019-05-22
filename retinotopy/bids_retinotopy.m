%% BIDS-based retinotopy solver for the project

base_path = '/Volumes/server/Projects/Retinotopy/CMag/data/BIDS';
cd(base_path);

vista_path = fullfile(base_path, 'derivatives', 'vistasoft');

addpath(genpath(fullfile(vista_path, 'shared', 'code')));

d = dir(base_path);
for ii = 1:numel(d)
    flnm = d(ii).name;
    if ~startsWith(flnm, 'sub-') || ~d(ii).isdir, continue; end
    sub = flnm(5:end);
    if ~strcmp(sub, 'wlsubj085'), continue; end
    fprintf('\n\n *** Processing subjects %s ***\n\n', sub);
    % these are the same for all subjects
    ses = 'nyu3t01'; 
    task = 'prf';
    % setup the subject!
    vista_path = fullfile('derivatives', 'vistasoft', ['sub-' sub], ['ses-' ses]);
    if ~exist(vista_path, 'dir')
        try
            bidsInitVista(base_path, sub, ses, task, 'preprocessed', vista_path);
        catch err
            warning('Error initializing vistasoft for %s / %s\n', sub, ses);
            rethrow(err);
            continue;
        end
    end
    cd(base_path);
    % now solve the prfs
    try
        solve_pRFs(base_path, sub, ses);
    catch err
        warning('Error solving pRFs for %s / %s\n', sub, ses);
    end
    fprintf('\n\n === Subjects %s complete ===\n\n', sub);
end


