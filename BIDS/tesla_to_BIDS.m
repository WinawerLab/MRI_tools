function tesla_to_BIDS(sourceDir, destDir, subj, varargin)
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
% sourceDir    = '/Volumes/server/Projects/BAIR/MRI/data/visual/wl_subj001_20170823/RAW';
% destDir      = tempdir;
% project      = 'visualFullSet';
% subj         = 'wl_subj001';
% session      = 'nyu3T01';
% distortScans = [5 6];
% distortDirs  = {'AP' 'PA'}
% epis         = 8:2:30;
% sbref        = 7:2:29;
% task         = {'hrf1' 'hrf2' 'hrf3' 'hrf4' 'hrf5' 'hrf6' 'task1' 'task2' ...
%                   'task3' 'task4' 'task5' 'task6'};
% tsvFiles     = {'sub-wlsubj001_ses-nyu3T01_task-hrf_run-01_bold.tsv' ...
%               'sub-wlsubj001_ses-nyu3T01_task-hrf_run-02_bold.tsv' ...
%               'sub-wlsubj001_ses-nyu3T01_task-hrf_run-03_bold.tsv' ...
%               'sub-wlsubj001_ses-nyu3T01_task-hrf_run-04_bold.tsv' ...
%               'sub-wlsubj001_ses-nyu3T01_task-hrf_run-05_bold.tsv' ...
%               'sub-wlsubj001_ses-nyu3T01_task-hrf_run-06_bold.tsv'};
% stimFile     = ...
% '/Users/jonathanwinawer/matlab/toolboxes/vistadispBAIR/Retinotopy/storedImagesMatrices/hrf_fMRI_1.mat';
%
% tesla_to_BIDS(sourceDir, destDir, subj, 'session', session, 'distortScans', distortDirs, ...
%   'epis', epis, 'sbref', sbref, 'task', task, 'tsvFiles', tsvFiles, ...
%   'stimFile', stimFile);
    
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
p.addRequired('project',@ischar);
p.addRequired('session', @(x) contains(x, 'nyu3T'));


p.addOptional('distortScans', @isnumeric); % vector of run numbers with distortion scans 
p.addOptional('distortDirs',   @iscell);   % cell array with PE directions for distortion scans

% For each EPI scan, we need several inputs: sbref, task name, tsvfile
p.addOptional('epis',  @isnumeric);    % vector of run numbers with EPI scans
p.addOptional('sbref', @isnumeric);    % vector of run numbers with sbref scans 
p.addOptional('task',  @iscell);       % cell array with one task name per EPI
p.addOptional('tsvFiles', @(x) iscell(x) || exist(x, 'file'));

p.addOptional('stimFile', @(x) iscell(x) || exist(x, 'file'));


p.parse(varargin{:});



end