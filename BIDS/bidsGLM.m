function results = bidsGLM(projectDir, subject, session, tasks, runnums, ...
        dataFolder, dataStr, designFolder, stimdur, modelType, glmOptsPath, tr)
%
% results = bidsGLM(projectDir, subject, [session], [tasks], [runnums], ...
%        [dataFolder], [dataStr], [designFolder],[stimdur], [modelType],[glmOptsPath], [tr]);
%
% Input
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
%     dataStr:          text string to specify filename for data 
%                           default: 'bold';
%     designFolder:     folder name containing design mats for glmSingle
%                           This specifies that the design matrices
%                           for each subject reside in
%                               <projectDir>/derivatives/design_matrices/<designFolder>/<subject>/<session>/
%                           and should contain either .tsv or .mat files
%                           default = [], which means no subfolder:
%                               <projectDir>/derivatives/design_matrices/<subject>/<session>/
%     stimdur:          duration of trials in seconds
%                           default = tr;
%     modelType:        name of folder to store outputs of GLMsingle (string)
%                           default = designFolder;
%     glmOptsPath:      path to json file specifying GLMsingle options
%                           default = [];
%     tr:               repetion time for EPI
%                           default = same as TR from raw data EPI
% Output
%     results:          structured array with GLMsingle results
%                           See GLMestimatesingletrial for details.
%
% Dependencies
%     GLMsingle repository (https://github.com/cvnlab/GLMsingle)
%     
%
% Example 1
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual_BIDS_compatible'; 
%     subject           = 'umcuchaam';
%     session           = 'umcu3tday139';
%     tasks             = 'temporalpattern';
%     runnums           = 1:4;
%     dataFolder        = 'fmriprep'; 
%     dataStr           = 'fsnative';
%     designFolder      = 'temporalpattern';
%     stimdur           = 0.5;
%     modelType         = 'temporalpatternRoundedTR';
%     glmOptsPath       = '~/matlab/toolboxes/BAIRanalysis/files/glmOpts.json';    
%     tr                = .85;
%       
%     % make the design matrices
%     bidsTSVtoDesign(projectDir, subject, session, tasks, runnums, designFolder);
%     % run the GLM
%     bidsGLM(projectDir, subject, session, tasks, runnums, ...
%         dataFolder, dataStr, designFolder, stimdur, modelType, glmOptsPath, tr)
%
% Example 2
%   Two step GLM: (1) optimze the hRF from an hrf-specific scan, then use the
%           optimzed hRF for denoising other scan types
% 
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%     subject           = 'wlsubj051';
%     session           = 'nyu3t01';
%     tasks             = 'hrfpattern';
%     runnums           = [];
%     dataFolder        = 'fmriprep';
%     dataStr           = 'fsnative';
%     designFolder      = 'hrfRoundedTROptimize';
%     stimdur           = 0.5;
%     modelType         = 'hrfRoundedTROptimize';
%     glmOptsPath       = 'glmOptsOptimize.json';
% 
%     % FIRST GLM
%     % make the design matrices
%     bidsTSVtoDesign(projectDir, subject, session, tasks, runnums, designFolder);
%     % run it
%     bidsGLM(projectDir, subject, session, tasks, runnums, ...
%         dataFolder, dataStr, designFolder, stimdur, modelType, glmOptsPath)
%     % look at the hrf
%     figure, plot(results.models{1})
% 
%     % SECOND GLM
%     tasks             = {'spatialobject' 'spatialpattern' 'temporalpattern'};
%     runnums           = [];
%     designFolder      = 'roundedTR';
%     glmOptsPath       = 'glmOptsAssume.json';
%     json = jsondecode(fileread(glmOptsPath));
%     json.hrfknobs     = results.models{1};
%     glmOptsPath       = fullfile(tempdir,glmOptsPath);
%     savejson('', json, 'FileName', glmOptsPath);
% 
%     bidsTSVtoDesign(projectDir, subject, session, tasks, runnums, designFolder);
%     % run it
%     bidsGLM(projectDir, subject, session, tasks, runnums, ...
%         dataFolder, dataStr, designFolder, stimdur, modelType, glmOptsPath)
% Example 3
%    GLM with temporally upsampled data/design matrices: 
% 
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%     subject           = 'wlsubj051';
%     session           = [];
%     tasks             = {'spatialobject' 'spatialpattern' 'temporalpattern'};
%     runnums           = [];
%     dataFolder        = 'threeTasksUpsampled';
%     designFolder      = 'threeTasksUpsampled';
%     stimdur           = [];
%     modelType         = [];
%     tr                = .85/5 ; % original data upsampled by 5x
%     glmOptsPath       = [];
%     dataStr           = 'fsnative';
%
%     % make the design matrices
%     bidsTSVtoDesign(projectDir, subject, session, tasks, runnums, designFolder, tr, dataFolder);
%
%     bidsGLM(projectDir, subject, session, tasks, runnums, ...
%       dataFolder, dataStr, designFolder, stimdur, modelType, glmOptsPath, tr)
%
%% Check inputs

if ~exist('session', 'var'),     session = [];      end
if ~exist('tasks', 'var'),       tasks   = [];      end
if ~exist('runnums', 'var'),     runnums  = [];     end
if ~exist('glmOptsPath', 'var'), glmOptsPath = [];  end
if ~exist('dataStr', 'var'),     dataStr = 'bold';  end

[session, tasks, runnums] = bidsSpecifyEPIs(projectDir, subject,...
    session, tasks, runnums);

% <dataFolder>
if ~exist('dataFolder', 'var') || isempty(dataFolder)
    dataFolder = 'preprocessed';
end
dataPath = fullfile (projectDir,'derivatives', dataFolder,...
    sprintf('sub-%s',subject), sprintf('ses-%s',session));
if exist(fullfile(dataPath, 'func'), 'dir')
    dataPath = fullfile (dataPath, 'func');
end
rawDataPath = fullfile(projectDir, sprintf('sub-%s', subject), ...
    sprintf('ses-%s', session), 'func');

% <designFolder>
if ~exist('designFolder', 'var'), designFolder = []; end
designPath = fullfile(projectDir, 'derivatives', 'design_matrices', ...
    designFolder, sprintf('sub-%s',subject), sprintf('ses-%s',session));
if ~exist(designPath, 'dir')
    warning('Subjectwise design path not found: %s', designPath); 
    
    disp('Checking for a design path shared across subjects..');
    designPath = fullfile(projectDir, 'derivatives', 'design_matrices', ...
    designFolder);
    if ~exist(designPath, 'dir')
        error('Shared design path not found: %s', designPath); 
    end
end

% <modelType>
if ~exist('modelType', 'var') || isempty(modelType)
    modelType = designFolder;
end

disp('design folder')
designPath

disp('Done checking inputs in bidsGLM')

%% Create GLMestimatesingletrial inputs

%****** Required inputs to GLMsingle *******************
% < design>
design = getDesign(designPath, tasks, runnums);

% <data>
data = bidsGetPreprocData(dataPath, dataStr, tasks, runnums);

% <tr>
if ~exist('tr', 'var') || isempty(tr)
    tr = bidsGetJSONval(rawDataPath,tasks, runnums, 'RepetitionTime');
    tr = cell2mat(tr);
    if length(unique(tr)) > 1
        disp(unique(tr))
        % TO DO: check if this is true for GLMsingle
        error(['More than one TR found:' ...
            'GLMdenoise expects all scans to have the same TR.'])
    else
        tr = unique(tr);
    end
end
% <stimdur>
if ~exist('stimdur', 'var') || isempty(stimdur), stimdur = tr;  end

%****** Optional inputs to GLMsingle *******************
% glm opts
opt = getGlmOpts(glmOptsPath);


%   <figuredir>
for ti=1:numel(tasks)
    figuredir   = fullfile (projectDir,'derivatives','GLMsingle', modelType, ...
                     sprintf('sub-%s',subject), sprintf('ses-%s',session), tasks{ti});

    if ~exist(figuredir, 'dir'); mkdir(figuredir); end
end


%% Save input arguments

inputVar = struct('projectDir', projectDir, 'subject', subject, ...
    'session', session, 'tasks', tasks, 'runnums', runnums, ...
    'dataFolder', dataFolder, 'dataStr', dataStr, 'designFolder', designFolder, ...
    'stimdur', stimdur, 'modelType', modelType, 'glmOptsPath', glmOptsPath, 'tr', tr);
    
fname = sprintf('sub-%s_ses-%s_%s_inputVar.json', subject, session, modelType);

savejson('',inputVar,fullfile(figuredir,fname));


%% Run the algorithm

% Cross validation criteria are checked within the function (e.g., checks
% if repeat trials occur across runs to do the cross validation)

[results,resultsdesign] = GLMestimatesingletrial(design,data,stimdur,tr,figuredir,opt);

% save the results
fname = sprintf('sub-%s_ses-%s_%s_results', subject, session, modelType);
save(fullfile(figuredir, fname), 'results', 'resultsdesign', '-v7.3');

end


%% ******************************
% ******** SUBROUTINES **********
% *******************************
function design = getDesign(designPath, tasks, runnums)
% <design> is the experimental design.  There are three possible cases:
%   1. A where A is a matrix with dimensions time x conditions.
%      Each column should be zeros except for ones indicating condition onsets.
%      (Fractional values in the design matrix are also allowed.)
%   2. {A1 A2 A3 ...} where each of the A's are like the previous case.
%      The different A's correspond to different runs, and different runs
%      can have different numbers of time points.
%   3. {{C1_1 C2_1 C3_1 ...} {C1_2 C2_2 C3_2 ...} ...} where Ca_b is a vector of
%      onset times (in seconds) for condition a in run b.  Time starts at 0
%      and is coincident with the acquisition of the first volume.  This case
%      is compatible only with <hrfmodel> set to 'assume'.
%   Because this function involves cross-validation across runs, there must
%   be at least two runs in <design>.
%

    tsvFiles = dir([designPath '/*_design.tsv']);
    matFiles = dir([designPath '/*_design.mat']);
    
    if ~isempty(tsvFiles) % first check for .tsv files
        fileNames = {tsvFiles.name}; listFiles = tsvFiles;
    elseif ~isempty(matFiles) % otherwise check for .mat files
        fileNames = {matFiles.name}; listFiles = matFiles;
    else
        error('No tsvs or mats found in: %s', designPath)
    end
    
    % this assumes each run is saved separately
    scan = 1; sharedDesignm = 0;
    for ii = 1:length(tasks)
        for jj = 1:length(runnums{ii})
            
            % check for both 0-padded and non-0-padded integers, throw
            % an error if we find other than 1.
            fileIdx = contains(fileNames,sprintf('task-%s_run-%d_',...
                tasks{ii}, runnums{ii}(jj)));
            fileIdx = or(fileIdx, contains(fileNames,sprintf('task-%s_run-%02d_',...
                tasks{ii}, runnums{ii}(jj))));
            if sum(fileIdx) > 1
               error(['Found more than one design matrix for task %s, run %d / '...
                      '%02d'], tasks{ii}, runnums{ii}(jj), runnums{ii}(jj))
            elseif sum(fileIdx) == 0 
                warning (['Subjectwise design Matrix *_task-%s_run-%d_design.tsv'...
                        ' was not found'],tasks{ii}, runnums{ii}(jj))
                disp('Will attempt to use a master file with design matrix:')
                sharedDesignm = 1;
                break;
            else
                % if tsv
                design{scan} = load(fullfile(listFiles(fileIdx).folder,listFiles(fileIdx).name));
                
                design_vars = load(fullfile(listFiles(fileIdx).folder,listFiles(fileIdx).name));
                design{scan} = design_vars.dms;

                scan         = scan+1;
            end
        end
    end

    % this assumes each run is saved as one aggregate file, where the
    % design matrix is a cell array, with each cell as a run.
    if sharedDesignm
        design = cell(1, length(tasks));
        for ii = 1:length(tasks)
            design_vars = load(fullfile(designPath, sprintf('%s_design', tasks{ii})));
            design_dict{ii} = design_vars.dms;

            % create a new design to match the included runs
            currentTaskdm = cell(1, length(runnums{ii}));
            
            % Assign values based on runnums{ii} indices
            for k = 1:length(runnums{ii})
                
                if any(runnums{ii}>length(design_dict{1}))
                    error(['Found run numbers that are outside of the index range  / '...
                        'of design matrix array %s'], fullfile(designPath, sprintf('%s_design', tasks{ii})))
                end

                currentTaskdm{k} = design_dict{1}{runnums{ii}(k)};
            end
            design{ii} = currentTaskdm;
        end
        design = [design{:}]; 
    end

end


function opt = getGlmOpts(glmOptsPath)

    if ~exist('glmOptsPath', 'var') || isempty(glmOptsPath)
        glmOptsPath = glmOptsMakeDefaultFile; 
    end
    
    json = jsondecode(fileread(glmOptsPath));
    
    if isfield(json, 'opt'), opt = json.opt; else, opt = []; end

end


function pth = glmOptsMakeDefaultFile()
    % see GLMestimatessingletrial for descriptions of optional input 

    json.opt.wantlibrary = 1;               % default=1 
    json.opt.wantglmdenoise = 1;            % default=1 
    json.opt.wantfracridge = 1;             % default=1 
    json.opt.chunknum = [];                 % default=50000
    json.opt.xvalscheme = [];               % default=leave one run out crossval
    json.opt.sessionindicator = [];         % default=ones (same session)

    % i/o flags

    json.opt.wantfileoutputs = [1,1,1,1];   % default=[1,1,1,1]
    json.opt.extraregressors = [];          % default=none
    json.opt.maxpolydeg = [];               % default=round(runMinutes/2) 
    json.opt.wantpercentbold = 1;           % default=1
    json.opt.hrftoassume = [];              % default=used only for ON/OFF
    json.opt.hrflibrary = [];               % default is matrix of 20 hRFs from getcanonicalhrflibrary.m
    json.opt.firdelay = [];                 % default=30
    json.opt.firpct = [];                   % default=99

    % model B flags

    json.opt.wantlss = [];                  % default=0 (uses ordinary least squares method)
    json.opt.numpcstotry = [];              % default=10 (maximum principal components)
    json.opt.brainthresh = [];              % default=[99 0.1]
    json.opt.brainR2 = [];                  % default is automatically determined
    json.opt.brainexclude = [];             % default=0 (any voxels can be included)
    json.opt.pcR2cutoff = [];               % default is automatically determined
    json.opt.pcR2cutoffmask = [];           % default=1 (any voxels can be included)   
    json.opt.pcstop = [];                   % default=1.05

    % model D flags

    json.opt.fracs = [];                    % default=fliplr(.05:.05:1)
    json.opt.wantautoscale = [];            % defaults=1

                  
    pth = fullfile(tempdir, 'glmOpts.json');
    savejson('', json, 'FileName', pth);
end          
           