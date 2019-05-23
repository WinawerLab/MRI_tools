function bids_export_niftis(projectPath, subject, session, task, freesurfer_id, fs_subjects_dir, output_path)
%% Configuration
    if startsWith(subject, 'sub-')
        subject = subject(5:end);
    end
    if startsWith(session, 'ses-')
        session = session(5:end);
    end
    if startsWith(task, 'task-')
        task = task(6:end);
    end
    % If these are set, then they won't be deduced from the script location in
    % the next section:
    % sub_name = 'wl_subj001';
    % session_path = '/Volumes/server/Projects/Retinotopy/wl_subj001/20171031_ColorRetinotopy';
    % session_name = '20171031_ColorRetinotopy';
    % stimuli_path = fullfile(session_path, 'Stimuli');
    % output_path = fullfile(session_path, 'Outputs');

    % If you don't have SUBJECTS_DIR set, then you'll need to set this manually
    % to be your FreeSurfer subject's directory
    % fs_subjects_dir = '/Volumes/server/Freesurfer_subjects';

    % If your subject has a different FreeSurfer subject ID than VistaSoft ID,
    % you must set this to the subject's freesurfer id:
    % freesurfer_id = 'wl_subj001';

    %% Deducing the remaining configuration data
    %  We assume that this file is in <something>/<subject>/code/scriptname.m
    %  for the purposes of deducing various paths.
    ssub = sprintf('sub-%s', subject);
    sses = sprintf('ses-%s', session);
    session_path = fullfile(projectPath, 'derivatives', 'vistasoft', ssub, sses)

    fprintf('Subject: %-12s  Session: %-20s\n', subject, session);

    % Next, figure out freesurfer data if not given
    if ~exist('freesurfer_id', 'var'), freesurfer_id = subject; end
    if ~exist('fs_subjects_dir', 'var')
        fs_subjects_dir = getenv('SUBJECTS_DIR');
    end

    % Last, figure out the output directory if not given
    if ~exist('output_path', 'var')
        output_path = fullfile(session_path, 'Outputs');
    end
    if ~exist(output_path, 'dir')
        sprintf('output_path (%s) not found, creating it', output_path);
        mkdir(output_path)
    end


    %% Navigate and initialize session

    cd(session_path);

    % The datsets we need to deal with:  (#output_plan)
    % This can be edited to ensure that various datasets are output with the
    % correct name. By default, exports the 'Full' datasets as a set of nifti
    % files with the prefix 'full'. To add more datasets, follow the example
    % of the scan_plan variable in the solve_pRFs script; here is an exmaple
    % for exporting a validation dataset and several training datasets:
    % output_plan.Test = 'test';
    % output_plan.Trn1 = 'trn1';
    % output_plan.Trn2 = 'trn2';
    % output_plan.Trn3 = 'trn3';
    % The above would output the datasets named 'Test', 'Trn1', 'Trn2', and
    % 'Trn3' as sets of files with the prefix 'test', 'trn1', 'trn2', and
    % 'trn3' respectively.
    output_plan = struct();
    output_plan = setfield(output_plan, task, 'full');

    % organize that...
    datasets = fieldnames(output_plan);
    outnames = cellfun(@(nm)(getfield(output_plan, nm)), datasets, 'UniformOutput', 0);

    % Do each in turn...
    vw = initHiddenGray();
    for ii = 1:numel(datasets)
        vw = viewSet(vw, 'current dt', datasets{ii});
        % load a model
        vw = rmSelect(vw, true, sprintf('Gray/%s/rm_%s-fFit.mat', datasets{ii}, datasets{ii}));
        rm = viewGet(vw, 'retinotopy model');

        % Load and export the angle/eccen
        vw = rmLoad(vw, 1, 'x0', 'map');
        vw = viewSet(vw, 'displaymode', 'map');
        functionals2nifti(vw, [], sprintf('%s/%s-xcrds.nii.gz', output_path, outnames{ii}));
        vw = rmLoad(vw, 1, 'y0', 'map');
        functionals2nifti(vw, [], sprintf('%s/%s-ycrds.nii.gz', output_path, outnames{ii}));
        % we also need the variance explained and sigma
        vw = rmLoad(vw, 1, 'sigma', 'map');
        functionals2nifti(vw, [], sprintf('%s/%s-sigma.nii.gz', output_path, outnames{ii}));
        vw = rmLoad(vw, 1, 'variance explained', 'map');
        functionals2nifti(vw, [], sprintf('%s/%s-vexpl.nii.gz', output_path, outnames{ii}));
    end

    % That's all!
end