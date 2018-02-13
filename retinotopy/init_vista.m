%% How this file works:
%  (1) You create a VistaSoft subject directory (but not an anatomy
%      directory)
%  (2) You process your subject with FreeSurfer
%  (3) You make Raw, Preproc, and Code directories in your VistaSoft dir
%  (4) You download your raw data to the Raw directory and run it with
%      Serra's preproc python script (such that the Preproc directory is
%      populated)
%  (5) You run this script to initialize anatomy and VistaSoft data

% Okay, init vistasoft
tbUse('vistasoft');

%% Configuration

% This directory should hold the path to the anatomies for individual
% subjects:
anatomies_path = '/Volumes/server/Projects/Anatomy';

% If these are set, then they won't be deduced from the script location in
% the next section:
% sub_name = 'wl_subj001';
% session_path = '/Volumes/server/Projects/Retinotopy/wl_subj001/20171031_ColorRetinotopy';
% session_name = '20171031_ColorRetinotopy';

% If you don't have SUBJECTS_DIR set, then you'll need to set this manually
% to be your FreeSurfer subject's directory
% fs_subjects_dir = '/Volumes/server/Freesurfer_subjects';

% If your subject has a different FreeSurfer subject ID than VistaSoft ID,
% you must set this to the subject's freesurfer id:
% freesurfer_id = 'wl_subj001';

% The annot_pattern is the string used to create the scan names via the
% sprintf function; e.g., if you have 4 scans that are all pRF scans, and
% you choose 'myPRF%02d' for this string,
% then your scans will be named myPRF01, myPRF02, myPRF03, and myPRF04.
annot_pattern = 'PRF%02d';

%% Deducing the remaining configuration data
%  We assume that this file is in <something>/<subject>/code/scriptname.m
%  for the purposes of deducing various paths.

% First, get the session path and subject name
if ~exist('sub_name', 'var')
    script_path = mfilename('fullpath');
    disp(script_path);
    [code_path,     script_name,  ~] = fileparts(script_path);
    [session_path,  code_name,    ~] = fileparts(code_path);
    [sub_path,      session_name, ~] = fileparts(session_path);
    [subjects_path, sub_name,     ~] = fileparts(sub_path);

    if ~strcmpi(code_name, 'code')
        error(['If initialization script is not in <subj>/<sess>/code'
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

%% Initializing the 3DAnatomy directory

% Whole brain reference
anatomy_path = fullfile(session_path, '3DAnatomy');
if exist(anatomy_path, 'dir')
    fprintf('Found anatomy directory %s...\n', anatomy_path);
else
    ext_anatomy_path = fullfile(anatomies_path, sub_name);
    if ~exist(ext_anatomy_path, 'dir')
        fprintf('Initializing anatomy data in %s...\n', ext_anatomy_path);
        init_anatomy(ext_anatomy_path, freesurfer_id, fs_subjects_dir);
    end
    fprintf('Linking anatomy directory %s to %s...\n', ...
           ext_anatomy_path, anatomy_path);
    cmd = sprintf('ln -s %s %s', ext_anatomy_path, anatomy_path);
    [sys_sts, sys_res] = system(cmd);
    if sys_sts ~= 0
        error(sprintf('Could not link anatomy dir %s to %s;\n%s\n', ...
                      ext_anatomy_path, anatomy_path, sys_res));
    end
end

t1_file = fullfile(anatomy_path, 't1.nii.gz');
if ~exist(t1_file,'file'), error('Anatomy file %s not found', t1_file); end


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
fprintf('Found %d preprocessed EPI files...\n', numel(epi_files));

% Inplane (mean of unwarped, merged, distortion scan)
ip = fullfile(preproc_path, 'distortion_merged_corrected_mean.nii.gz');
if ~exist(ip, 'file'), error('Inplane file %s not found', ip); end

% alignment file
align_file = fullfile(preproc_path, 'distort2anat_tkreg.dat');
if ~exist(align_file, 'file')
    error(sprintf('Alignment file %s not found', align_file));
end



%% Generating the params structure

% trim off trailing / of the session-path if need-be:
while session_path(end) == filesep, session_path = session_path(1:end-1); end

% move ourselves over to the session's directory
cd(session_path);

% remove absolute paths where appropriate so that the pRF solutions can be run
% on the HPC
spl = numel(session_path);
if startsWith(ip, session_path)
    ip = ip(spl+2:end);
end
for ii = 1:numel(epi_files)
    if startsWith(epi_files{ii}, session_path)
        epi_files{ii} = epi_files{ii}(spl+2:end);
    end
end
if startsWith(t1_file, session_path)
    t1_file = t1_file(spl+2:end);
end

disp(ip);
disp(t1_file);

params = mrInitDefaultParams;

% And insert the required parameters: 
params.inplane      = ip;
params.functionals  = epi_files;
params.sessionDir   = session_path;

% Specify some optional parameters
params.vAnatomy = t1_file;
params.subject  = sub_name;

% Name each scan
params.annotations = cell(numel(epi_files), 1);
for ii = 1:numel(epi_files)
    params.annotations{ii} = sprintf(annot_pattern, ii);
end

%% Running the initialization 
ok = mrInit(params); 

% open a session without GUI
userpathvw = initHiddenInplane();

fid = fopen(align_file, 'r');

% skip header lines
for k=1:4, tline = fgets(fid); end
for k = 1:4; R(k,:) = str2num(fgetl(fid)); end
fclose(fid);
vistaAlignment = fs_tkreg2vista(R, ip, t1_file);
mrSESSION = sessionSet(mrSESSION, 'alignment', vistaAlignment);
saveSession();

% Install gray segmentation
% path to class file
t1class_path = fullfile(anatomy_path, 't1_class.nii.gz');
% number of layers in gray graph along surface (3 layers for 1 mm voxels)
numGrayLayers = 3;
installSegmentation([], [], t1class_path, numGrayLayers);

% Import freesurfer meshes
lsrf_file = fullfile(anatomy_path, 'Left', '3DMeshes', 'Left_white.mat');
if ~exist(lsrf_file, 'file')
    fprintf('Initializing FreeSurfer meshes...\n');
    sub_fs_path = fullfile(fs_subjects_dir, freesurfer_id);
    meshImportFreesurferSurfaces(sub_fs_path, 'b');
end

%% Helper Functions

% Helper function for initializing anatomy from FreeSurfer
function init_anatomy(anat_path, fs_subj_id, fs_subjects_dir)
    % see if we have a subjects dir
    if isempty(fs_subjects_dir)
        error('no SUBJECTS_DIR given');
    elseif ~exist(fs_subjects_dir, 'dir')
        error(sprintf('SUBJECTS_DIR not found: %s', fs_subjects_dir));
    end
    data_dir = fullfile(fs_subjects_dir, fs_subj_id);
    if ~exist(data_dir, 'dir')
        error(sprintf('FreeSurfer directory (%s) not found', data_dir));
    end
    % Make the VistaSoft anatomy directory and check it
    mkdir(anat_path);
    if ~exist(anat_path, 'dir')
        error(sprintf('could not make anatomy directory %s', anat_path));
    end
    outfile = fullfile(anat_path, 't1_class.nii.gz');
    t1file = fullfile(anat_path,  't1.nii.gz');
    alignTo = fullfile(data_dir,  'mri', 'orig.mgz');
    fs_ribbon2itk(fs_subj_id, outfile, true, alignTo);
    % Check that you created a t1 class file (ribbon) and t1 anatomy
    if ~exist(outfile, 'file'),
        error(sprintf('Failed to create class file %s', outfile));
    end
    if ~exist(t1file, 'file'),
        error(sprintf('Failed to create T1 file %s', t1file));
    end
end