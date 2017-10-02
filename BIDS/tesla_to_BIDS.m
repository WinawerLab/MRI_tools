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

error('Function not yet finished')
    
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
p.addRequired('datadir',exist(x, 'dir'));  % source directory
p.addRequired('outdir');                   % destination directory
p.addRequired('project', 'visualFullSet',@ischar);
p.addRequired('subject', 'wlsubj000' ,@ischar);
p.addRequired('session', 'nyu3T01' ,@(x) contains(x, 'nyu3T'));


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