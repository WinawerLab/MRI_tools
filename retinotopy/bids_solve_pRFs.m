function bids_solve_pRFs(projectPath, subject, session, task, fsid, fssubdir, stim_dir)
% solve_pRFs() averages the prf scans (with the given taskname, assumed to
% be 'prf') for the given subject and session in the given BIDS project
% directory.

path0 = pwd();
ssub = sprintf('sub-%s', subject);
sses = sprintf('ses-%s', session);
cd(projectPath);
try
    wd = pwd();
    
    dataset_name = task;
    session_path = fullfile('derivatives', 'vistasoft', ssub, sses);
    output_path  = fullfile(session_path, 'Gray', dataset_name);
    assert(exist(session_path, 'dir') ~= 0);

    if ~exist('fsid', 'var') || isempty(fsid), fsid = sub; end
    assert(~isempty(fsid));
    if ~exist('fssubdir', 'var') || isempty(fssubdir)
        fssubdir = getenv('SUBJECTS_DIR');
    end
    if ~exist('stim_dir', 'var') || isempty(stim_dir)
        stim_dir = fullfile(path0, 'files');
    end
    assert(~isempty(fssubdir));
    assert(~isempty(fullfile(fssubdir, fsid)));

    fprintf('Subject: %-12s  Session: %-20s\n', subject, session);

    % stim params and such should be in MRI_tools/retinotopy/files directory
    params_flnm = fullfile(stim_dir, 'scan_params.mat');
    images_flnm = fullfile(stim_dir, 'scan_images.mat');
    if ~exist(params_flnm, 'file') || ~exist(images_flnm, 'file')
        error('Could not file params or images file');
    end

    %% Navigate and initialize session
    cd(session_path);
    mrvCleanWorkspace();

    gr = initHiddenGray();
    vw = initHiddenInplane();
    global dataTYPES;

    % First off, if the model file already exists, skip this all...
    ffit_path = fullfile(wd, output_path);
    ffit_flnm = fullfile(ffit_path, sprintf('rm_%s-fFit.mat', dataset_name));
    if ~exist(ffit_flnm, 'file')
        
        num_dts = numel(dataTYPES);

        scan_list = numel(dataTYPES(vw.curDataType).scanParams);

        vw = averageTSeries(vw, scan_list, dataset_name);

        ii = num_dts + 1;
        vw = viewSet(vw, 'current dt', ii);
        gr = viewSet(gr, 'current dt', ii);
        gr = ip2volTSeries(vw, gr, 0, 'linear');

        %% Initialize Retinotopy Parameters
        prfModels = 'one gaussian';
        load(params_flnm, 'params');

        % Set default retinotopy stimulus model parameters
        sParams = rmCreateStim(vw, gr);
        sParams.stimType   = 'StimFromScan'; % This means the stimulus images will
                                             % be read from a file.
        sParams.stimSize   = 12.4;           % stimulus radius (deg visual angle)
        sParams.nDCT       = 1;              % detrending frequeny maximum (cycles
                                             % per scan): 1 means 3 detrending
                                             % terms, DC (0 cps), 0.5, and 1 cps
        sParams.imFile     = images_flnm;    % file containing stimulus images
        sParams.paramsFile = params_flnm;    % file containing stimulus parameters
                                             % 'thresholdedBinary', whenreading in images, treat any
                                             % pixel value different from background as a 1, else 0
        sParams.imFilter   = 'thresholdedBinary';
        % we switch from the default positive Boynton hRF to the biphasic SPM style
        sParams.hrfType    = 'two gammas (SPM style)';
        % pre-scan duration will be stored in frames for the rm, but was stored in
        % seconds in the stimulus file
        sParams.prescanDuration = params.prescanDuration/params.framePeriod;


        %% Solve the retinotopy models
        gr = initHiddenGray();
        gr = viewSet(gr, 'current dt', ii);
        dataTYPES(ii).retinotopyModelParams = [];
        dataTYPES(ii) = dtSet(dataTYPES(ii), 'rm stim params', sParams);
        gr = rmMain(gr, [], 'coarse to fine', ...
                    'model', {prfModels}, ...
                    'matFileName', sprintf('rm_%s', dataTYPES(ii).name));
        % save the results
        saveSession();
    end

    vw = initHiddenGray();
    %% export the nifti files
    flnms = {'xcrds.nii.gz', 'ycrds.nii.gz', 'sigma.nii.gz', 'vexpl.nii.gz'};
    vsnms = {'x0', 'y0', 'sigma', 'variance explained'};
    vw = viewSet(vw, 'current dt', dataset_name);
    % load a model
    vw = rmSelect(vw, true, ffit_flnm);
    rm = viewGet(vw, 'retinotopy model');
    % Load and export the parameters
    for ii = 1:numel(flnms)
        disp(flnms{ii});
        flnm = fullfile(ffit_path, flnms{ii});
        if exist(flnm, 'file')
            fprintf('Skipping export of existing file %s\n', flnm);
            continue;
        else
            fprintf('Exporting file %s\n', flnms{ii});
            vw = rmLoad(vw, 1, vsnms{ii}, 'map');
            vw = viewSet(vw, 'displaymode', 'map');
            functionals2nifti(vw, [], flnm);
        end
    end    
    % return to original path
    cd(path0);
catch err
    cd(path0);
    rethrow(err)
end

