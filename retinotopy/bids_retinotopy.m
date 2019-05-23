function bids_retinotopy(bids_path, subject, session, task, anatomyDir)
% BIDS-based retinotopy solver

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
    % THIS IS NO GOOD
    % if ~exist(vista_path, 'dir')
        try
            bidsInitVista(bids_path, subject, session, task, 'preprocessed', vista_path, anatomyDir);
        catch err
            warning('Error initializing vistasoft for %s / %s\n', subject, session);
            rethrow(err);
        end
    % end
    cd(bids_path);
    % now solve the prfs
    % try
        bids_solve_pRFs(bids_path, subject, session);
    % catch err
    %     warning('Error solving pRFs for %s / %s\n', subject, session);
    % end
    fprintf('\n\n === Subject %s complete ===\n\n', subject);
end