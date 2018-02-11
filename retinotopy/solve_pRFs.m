%% Configuration

clear('sub_name', 'session_path', 'session_name', 'stimuli_path', ...
      'params_file', 'fs_subjects_dir', 'freesurfer_id');

% If these are set, then they won't be deduced from the script location in
% the next section:
% sub_name = 'wl_subj001';
% session_path = '/Volumes/server/Projects/Retinotopy/wl_subj001/20171031_ColorRetinotopy';
% session_name = '20171031_ColorRetinotopy';
% stimuli_path = fullfile(session_path, 'Stimuli');

% If you don't have SUBJECTS_DIR set, then you'll need to set this manually
% to be your FreeSurfer subject's directory
% fs_subjects_dir = '/Volumes/server/Freesurfer_subjects';

% If your subject has a different FreeSurfer subject ID than VistaSoft ID,
% you must set this to the subject's freesurfer id:
% freesurfer_id = 'wl_subj001';

%% Deducing the remaining configuration data
%  We assume that this file is in <something>/<subject>/code/scriptname.m
%  for the purposes of deducing various paths.

% First, get the session path and subject name
if ~exist('sub_name', 'var')
    script_path = mfilename('fullpath');
    [code_path,     script_name,  ~] = fileparts(script_path);
    [session_path,  code_name,    ~] = fileparts(code_path);
    [sub_path,      session_name, ~] = fileparts(session_path);
    [subjects_path, sub_name,     ~] = fileparts(sub_path);

    if ~strcmpi(code_name, 'code')
        error(['If initialization script is not in <subj>/<sess>/code' ...
               ' then it must be edited manually to include paths']);
    end
end
if ~exist('session_name', 'var') || ~exist('session_path', 'var')
    error('No session name/path deduced or provided');
end
fprintf('Subject: %-12s  Session: %-20s\n', sub_name, session_name);
   
% Next, figure out freesurfer data if not given
if ~exist('freesurfer_id', 'var'), freesurfer_id = sub_name; end
if ~exist('fs_subjects_dir', 'var')
    fs_subjects_dir = getenv('SUBJECTS_DIR');
end

% Figure out the stimulus path if not given
if ~exist('params_file', 'var') || ~exist('images_file', 'var')
    need_stimdir = true;
else
    need_stimdir = false;
end
if ~exist('stimuli_path', 'dir') && need_stimdir
    stimuli_paths = {'Stimulus', 'stimulus', 'Stimuli', 'stimuli', ...
                     'Stim', 'stim'};
    for i = 1:numel(stimuli_paths)
        stimuli_path = fullfile(session_path, stimuli_paths{i});
        if isdir(stimuli_path), break;
        else stimuli_path = [];
        end
    end
    if isempty(stimuli_path)
        error('Could not deduce the stimuli path')
    end
end
% And figure out the actual param and image files...
if ~exist('params_file', 'var')
    fls = dir(fullfile(stimuli_path, '*_params.mat'));
    params_file = [];
    for ii = 1:numel(fls)
        fldat = fls(ii);
        fl = fullfile(stimuli_path, fldat.name);
        if fl(1) == '.',  continue;
        elseif isdir(fl), continue;
        else              params_file = fl;
                          break;
        end
    end
    if isempty(params_file)
        error('Could not deduce location of parameters mat file');
    end
end
if ~exist(params_file, 'file')
    error(sprintf('params file (%s) does not exist', params_file));
end
if ~exist('images_file', 'var')
    fls = dir(fullfile(stimuli_path, '*_images.mat'));
    images_file = [];
    for ii = 1:numel(fls)
        fldat = fls(ii);
        fl = fullfile(stimuli_path, fldat.name);
        if fl(1) == '.',  continue;
        elseif isdir(fl), continue;
        else              images_file = fl;
                          break;
        end
    end
    if isempty(images_file)
        error('Could not deduce location of images mat file');
    end
end
if ~exist(images_file, 'file')
    error(sprintf('images file (%s) does not exist', images_file));
end

%% Finding the EPI files

% Find the preproc path...
preproc_paths = {'preproc', 'Preproc', 'preprocessed', 'Preprocessed'};
for i = 1:numel(preproc_paths)
    preproc_path = fullfile(session_path, preproc_paths{i});
    if isdir(preproc_path), break;
    else                    preproc_path = [];
    end
end
if isempty(preproc_path), error('Could not find a preproc directory');
else fprintf('Using preproc directory %s...\n', preproc_path);
end
d = dir(fullfile(preproc_path, 'timeseries_corrected*.nii.gz'));
for ii = 1:length(d)
    epi_files{ii} = fullfile(preproc_path, d(ii).name);
end

%% Initialize Retinotopy Parameters
prfModels = 'one gaussian';
load(params_file, 'params');

%% Navigate and initialize session

cd(session_path);

gr = initHiddenGray();
vw = initHiddenInplane();

num_dts = numel(dataTYPES);
num_scans = numel(epi_files);

% Average the retinotopy scans (#scan_plan)
% Edit this cell array to analyze alternate subsets of the scans;
% e.g., if you have 6 scans and want a few training and validation
% datasets, you might do something like:
% scan_plan.Full = 1:6 % the Full dataset
% scan_plan.Test = 1:3 % Test/Validation dataset
% scan_plan.Trn1 = [4] % Training dataset 1
% scan_plan.Trn2 = [5] % Training dataset 2
% scan_plan.Trn3 = [6] % Training dataset 3
scan_plan.Full = 1:num_scans; % the 'Full' dataset of all retinotopy scans

% iterate through these datasets and solve them
dsets = fieldnames(scan_plan);
for dset = 1:numel(dsets)
    name = dsets{dset};
    scan_list = getfield(scan_plan, name);
    
    % average it
    vw = averageTSeries(vw, scan_list, 'Full');

    ii = num_dts + dset;
    vw = viewSet(vw, 'current dt', ii);
    gr = viewSet(gr, 'current dt', ii);
    gr = ip2volTSeries(vw, gr, 0, 'linear');

    % Set default retinotopy stimulus model parameters
    sParams = rmCreateStim(vw, gr);
    sParams.stimType   = 'StimFromScan'; % This means the stimulus images will
                                         % be read from a file.
    sParams.stimSize   = 14;             % stimulus radius (deg visual angle)
    sParams.nDCT       = 1;              % detrending frequeny maximum (cycles
                                         % per scan): 1 means 3 detrending
                                         % terms, DC (0 cps), 0.5, and 1 cps
    sParams.imFile     = images_file;    % file containing stimulus images
    sParams.paramsFile = params_file;    % file containing stimulus parameters
    % 'thresholdedBinary',  whenreading in images, treat any pixel value
    % different from background as a 1, else 0
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
end

saveSession();

