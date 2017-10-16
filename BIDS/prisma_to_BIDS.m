function prisma_to_BIDS(sourceDir, destDir, subj, varargin)
% Some things we'll need to define in the outputs:
%   - Project name
%   - subject
%   - session name
%   - stimulius file (file in project/stimuli)
%   - tsv files for stimuli
%   - 3D Anatomy (if requested)
%   - freesurfer directory (if requested)
%   - task name for each EPI
%
% Example
%  sourceDir    = '/Volumes/server/Projects/BAIR/MRI/data/visual/wl_subj043_20170824/RAW';
%  destDir      = fullfile(tempdir, 'BIDSCompatibleOutput');
%  project      = 'visualFullSet';
%  subj         = 'wl_subj043';
%  session      = 'nyu3T01';
%  anat         = 5;
%  distortScans = [6 7];
%  distortDirs  = {'AP' 'PA'}
%  epis         = 9:2:29;
%  sbref        = 8:2:28;
%  task         = repmat({'spat'}, 1,11);
%  for ii = 1:11, tsvFiles{ii} = sprintf('sub-wlsubj043_ses-nyu3T01_task-spat_run-%02d_bold.tsv', ii); end        
%  stimFile     = '~/matlab/toolboxes/vistadispBAIR/Retinotopy/storedImagesMatrices/spatiotemporal_fMRI_1.mat';
%
%  prisma_to_BIDS(sourceDir, destDir, subj, 'session', session, 'distortScans', distortDirs, ...
%    'epis', epis, 'sbref', sbref, 'task', task, 'tsvFiles', tsvFiles, ...
%    'stimFile', stimFile, 'project', project, 'anat', anat);
    
%% Parse inputs
p = inputParser;

% Required:
%  datadir: Source directory
%  outdir:  Destination directory
%  project: Name of project (this is a string that is one level above the
%               subject)
%  subject: alphanumeric code for subject - could agree with freesurfer
%               directory name, as in wl_subj004, but will have to have
%               underscore removed for naming BIDS directories

% p.addRequired('sourceDir', @(x) exist(x, 'dir'));  % source directory
% p.addRequired('destDir');                   % destination directory
% p.addRequired('subj',@ischar);

% Project
p.addRequired('project',@ischar);

% Session 
p.addRequired('session', @(x) contains(x, 'nyu3T'));

% Field maps or distort scans
p.addOptional('distortScans', @isnumeric); % vector of run numbers with distortion scans 
p.addOptional('distortDirs',   @iscell);   % cell array with PE directions for distortion scans

% EPIS and related. For each EPI scan, we need several inputs: sbref, task name, tsvfile
p.addOptional('epis',  @isnumeric);    % vector of run numbers with EPI scans
p.addOptional('sbref', @isnumeric);    % vector of run numbers with sbref scans 
p.addOptional('task',  @iscell);       % cell array with one task name per EPI
p.addOptional('tsvFiles', @(x) iscell(x) || exist(x, 'file'));

% Stim file (path to stimulus file)
p.addOptional('stimFile', @(x) iscell(x) || exist(x, 'file'));

% Anatomy - either a scan number (if acquired from this session) or a freesurfer subjid
p.addOptional('anat', @(x) isnumeric(x) || ischar(x));

p.parse(varargin{:});

% folder Structures
projectDir = fullfile(destDir, p.Results.project);
subjectDir = fullfile(projectDir, sprintf('sub-%s', strrep(subj, '_', '')));
sessionDir = fullfile(subjectDir, sprintf('ses-%s', strrep(p.Results.session, '_', '')));


%% Go!
% Make folders
bidsCreateFolders(projectDir, subjectDir, sessionDir);

% Make anat folder and files
bidsMakeAnat(sourceDir, sessionDir, p);

% Make derivatives folder and files
% Make fmap folder and files
% Make func folder and files
% Make sourcedata folder and files

% Add stimulus files

% Add participant to partipicants.tsv

% Check BIDS compatibility


end

function bidsCreateFolders(projectDir, subjectDir, sessionDir)
% Check whether current project exists
if exist(projectDir, 'dir')
    % do nothing
else
    fprintf('Project %s does not yet exist. New project directory being created.\n', ...
        projectDir);
end


% Check whether current subject exists within this project
if exist(subjectDir, 'dir')
 % do nothing
else
    fprintf('Subject %s does not yet exist within this project.\n', subjectDir);
    fprintf('New subject directory being created.\n')
end

% Check whether current session exists within this project and session
if exist(sessionDir, 'dir')
    error('The folder %s already exists for this project, subject, and session. You must either delete the existing session directory or change the name of the input session directory.', ... 
        sessionDir);
else
    fprintf('New session directory %s being created.\n', sessionDir);
end

 
% parent folder
mkdir(sessionDir); 
end

function bidsMakeAnat(sourceDir, sessionDir, p)

anat = p.Results.anat;

anatDest = sprintf('sub-%s_ses-%s_T1w.%s', subj, session, extention);
fullfile(sessionDir, 'anat', 

if isnumeric(anat)
    % if it's a number, then it was an acquisition from this scan session
    anatFolder = dir(fullfile(sourceDir, sprintf('%02d*', anat)));
    switch numel(anatFolder)
        case 0
            error('Anatomy folder %s not found', fullfile(sourceDir, sprintf('%02d*', anat)))
        case 1
            anatFile = dir(fullfile(sourceDir, anatFolder.name, '*.nii*'));
        otherwise
            error('More than one anatomy folder found')
    end
        
    switch numel(anatFile)
        case 0
            error('Anatomy file not found in %s', fullfile(sourceDir, anatFolder.name))
        case 1
            anatSource = fullfile(sourceDir, anatFolder.name, anatFile.name);
        otherwise
            error('More than one anatomy file found in %s', fullfile(sourceDir, anatFolder.name))
    end
          
else
   % if it's a str, then assume it's a freesufer dirctory and find the
   % anatomy file
end


end