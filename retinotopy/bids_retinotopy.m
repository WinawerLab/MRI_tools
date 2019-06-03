function bids_retinotopy(bids_path, subject, session, task, anatomyDir, stim_dir)
% BIDS-based retinotopy solver

    if ~exist('anatomyDir', 'var') || isempty(anatomyDir)
        anatomyDir = fullfile('/Volumes', 'server', 'Projects', 'Anatomy');
    end
    if ~exist('stim_dir', 'var') || isempty(stim_dir)
        stim_dir = fullfile(pwd(), 'files');
    end

    
    if startsWith(subject, 'sub-')
        subject = subject(5:end);
    end
    if startsWith(session, 'ses-')
        session = session(5:end);
    end
    if startsWith(task, 'task-')
        task = task(6:end);
    end

    fprintf('\n\n *** Processing subject %s, session %s ***\n\n', subject, session);

    % setup the subject!
    vista_path = fullfile(bids_path, 'derivatives', 'vistasoft', ['sub-' subject], ['ses-' session]);
    if ~exist(fullfile(vista_path, 'mrInit_params.mat'), 'file')
        bidsInitVista(bids_path, subject, session, task, 'preprocessed', vista_path, anatomyDir);
    end
    % now solve the prfs
    bids_solve_pRFs(bids_path, subject, session, task, [], [], stim_dir);
    fprintf('\n\n === Subject %s complete ===\n\n', subject);
end