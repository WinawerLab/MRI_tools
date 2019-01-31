function bidsInitVista(projectDir, subject, session, tasks,... % runnums,
    dataFolder , analysisDir)
% function bidsInitVista(projectDir, subject, [session], [tasks], [runnums], ...
%    [dataFolder], [analysisDir])
%
%   Input:
%
%   Required
%
%     projectDir:       path where the BIDS projects lies (string)
%     subject:          BIDS subject name (string, all lower case)
%
%   Optional
%
%     session:          BIDS session name (string, all lower case)
%                           default: folder name inside subject dir (if
%                           more than one or less than one, return error)
%     tasks:            one or more BIDS tasks (string or cell array of strings)
%                           default: all tasks in session
%     runnums:           BIDS run numbers (vector or cell array of vectors)
%                           default: all runs for specified tasks
%     dataFolder:       folder name containing preprocessed BOLD data. Note
%                           that this folder should reside in
%                               <projectDir>/derivatives/
%                           and should contain subdirectories and BOLD
%                           files of the form
%                               <subject>/<session>/*.nii
%                           default: 'preprocessed'
%     analysisDir:      path to project directory containing non-BIDS analysis output
%                           default: Using projectDir, a project will be inferred
%                            fullfile(fileparts(projectDir), 'Analyses', <subj>, <session>);
% Example:
%
% projectDir  = '/Volumes/server/Projects/BAIR/Data/BIDS/tactilepilot';
% subject     = 'wlsubj063';
% analysisDir = '/Volumes/server/Projects/BAIR/Analyses/tactilepilot/sub-wlsubj063/ses-nyu3t01'
%
% bidsInitVista(projectDir, subject, [] , [], [], analysisDir)

%% Check inputs

if ~exist('session', 'var'), session = []; end
if ~exist('tasks', 'var'),   tasks = []; end 

[session, tasks, runnum] = bidsSpecifyEPIs(projectDir, subject, session, tasks); 

% <dataFolder>
if ~exist('dataFolder', 'var') || isempty(dataFolder)
    dataFolder = 'preprocessed';
end
dataPath = fullfile (projectDir,'derivatives', dataFolder,...
    sprintf('sub-%s',subject), sprintf('ses-%s',session));
assert(boolean(exist(dataPath, 'dir')))

% <Analysis folder>
if ~exist(analysisDir , 'dir')
    analysisDir = fullfile (fileparts (projectDir), 'Analyses', ...
        sprintf('sub-%s',subject), sprintf('ses-%s',session));
   
end
assert(boolean(exist(analysisDir, 'dir')))
cd(analysisDir);

%% Whole brain reference
% Set up file names and associated locations
anatFile  = 't1.nii.gz';
anatPath   = fullfile('3DAnatomy', anatFile);
classFile = 't1_class.nii.gz';
classPath  = fullfile('3DAnatomy', classFile);

% Look for anatomy files and make them if they don't exist
if ~exist(fullfile(analysisDir,'3DAnatomy'),'dir')
    anatomyPath = fullfile('/Volumes/server/Projects/Anatomy',subject);
    if ~exist(anatomyPath, 'dir') ||~exist(fullfile(anatomyPath, anatFile), 'file')...
            && ~exist(fullfile(anatomyPath,classFile), 'file')
        mkdir (anatomyPath);
        % make anatfile
        fshome = getenv('SUBJECTS_DIR');
        str = sprintf('mri_convert -i %s -o %s', fullfile (fshome, subject, 'mri','T1.mgz')...
            , fullfile(anatomyPath, anatFile));
        system(str)
        % make t1_classfile
        fs_ribbon2itk(subject, fullfile(anatomyPath, classFile));
    end
    % make sym link
    str = sprintf('ln -s %s %s', anatomyPath , '3DAnatomy');
    system(str)
end


%% Set params
cd (analysisDir);

% Generate the expected generic params structure
params = mrInitDefaultParams;

% Find and name each scan
d = dir(fullfile(dataPath, '*_preproc.nii.gz'));

% Need to add runnums and tasks
epiFiles = cell(1, length(d));

for ii = 1:length(d)
    epiFiles{ii} = d(ii).name;
    task   = bidsGet(d(ii).name,'task');
    runnum = bidsGet(d(ii).name,'run');
    params.annotations{ii} = sprintf('task-%s_run-%s', task, runnum);
end

% Inplane (mean of unwarped, merged, distortion scan)
ip = 'distortion_merged_corrected_mean.nii.gz';
if ~exist(fullfile(dataPath, ip), 'file')
    error('Inplane file %s not found.', fullfile(sessionDir, ip));
end

% And insert the required parameters:
params.inplane      = fullfile(dataPath, ip);
params.functionals  = fullfile(dataPath, epiFiles);
params.sessionDir   = analysisDir;

% Specify some optional parameters
params.vAnatomy = anatPath;
params.subject  = subject;

%% Run it:
ok = mrInit(params);

%% open a session without GUI
vw = initHiddenInplane();
mrGlobals();

% Install alignment
alignfile = fullfile(dataPath, 'distort2anat_tkreg.dat');
if ~exist(alignfile, 'file')
    error('Alignment file %s not found', alignfile)
end
fid = fopen(alignfile, 'r');

% skip header lines
for k=1:4, tline = fgets(fid); end
for k = 1:4; R(k,:) = str2num(fgetl(fid)); end
fclose(fid);
vistaAlignment = fs_tkreg2vista(R,params.inplane, fullfile('3DAnatomy',anatFile));
mrSESSION = sessionSet(mrSESSION, 'alignment', vistaAlignment);
saveSession();

% Install gray segmentation

% number of layers in gray graph along surface (3 layers for 1 mm voxels)
numGrayLayers = 3;
installSegmentation([], [], classPath , numGrayLayers);

% Import freesurfer meshes
if ~exist(fullfile('3DAnatomy', 'Left', '3DMeshes'), 'dir') || ...
        ~exist(fullfile('3DAnatomy', 'Right', '3DMeshes'), 'dir')
    meshImportFreesurferSurfaces(subject);
end

end

%% DEBUG
% % Check that initialization worked correctly
% 
% % Open an Inplane View
% vw = mrVista;
% 
% % Compute a mean map for the first scan
% vw=computeMeanMap(vw,getCurScan(vw));
% vw=loadMeanMap(vw, 1);
% vw=setDisplayMode(vw,'map');
% vw=refreshScreen(vw);
% vw=viewSet(vw, 'mapwin', [0.2 1]);
% 
% % If the mean functional data is aligned to the underlay, then we have good
% % agreement between EPIs and INPLANE
% 
% % Now we check that the gray matter is aligned to the EPIs
% vw=makeGrayROI(vw);
% vw=setDisplayMode(vw,'anat');
% vw=viewSet(vw, 'ROI color', 'k');
% vw=refreshScreen(vw);
% 
% % If you are not satisfied with the alignment, consider using the vistasoft
% % script s_alignInplaneToAnatomical to refine the alignment
% % s_alignInplaneToAnatomical;