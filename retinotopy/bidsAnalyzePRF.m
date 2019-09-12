function results = bidsAnalyzePRF(projectDir, subject, session, tasks, runnums, ...
        dataFolder, dataStr, apertureFolder, prfOptsPath, tr)
%
% results = bidsAnalyzePRF(projectDir, subject, [session], [tasks], [runnums], ...
%        [dataFolder], [apertureFolder], [glmOptsPath], [tr]);
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
%                           default: 'fmriprep'
%     dataStr:          text string to specify filename for data 
%                           default: 'bold';
%     apertureFolder:   folder name containing stimulus apertureas for analyzePRF
%                           This specifies that the stimulus apertures
%                           for each subject reside in
%                               <projectDir>/derivatives/stim_apertures/<apertureFolder>/<subject>/<session>/
%                           and should contain .mat files
%                           default = [], which means no subfolder:
%                               <projectDir>/derivatives/stim_apertures/<subject>/<session>/
%     prfOptsPath:      path to json file specifying analyzePRF options
%                           default = [];
%     tr:               repetion time for EPI
%                           default = same as TR from raw data EPI
% Output
%     results:          structured array with analyzePRF results
%                           See analyzePRF for details.
%
% Dependencies
%     analyzePRF repository (https://github.com/kendrickkay/analyzePRF)
%     
%
% Example 1
%     projectDir        = '/Volumes/server/Projects/SampleData/BIDS'; 
%     subject           = 'wlsubj042';
%     session           = '01';
%     tasks             = 'prf';
%     runnums           = 1:2;
%     dataFolder        = 'fmriprep'; 
%     dataStr           = 'fsnative*.mgz';
%     apertureFolder    = [];
%     prfOptsPath       = [];    
%     tr                = [];
%       
%     % make the stimulus apertures 
%     bidsStimulustoApertures(projectDir, subject, session, tasks, runnums, apertureFolder);
%
%     % run the prf analysis
%     bidsAnalyzePRF(projectDir, subject, session, tasks, runnums, ...
%        dataFolder, dataStr, apertureFolder, prfOptsPath, tr)
%

%% Check inputs

if ~exist('session', 'var'),     session = [];      end
if ~exist('tasks', 'var'),       tasks   = [];      end
if ~exist('runnums', 'var'),     runnums  = [];     end
if ~exist('prfOptsPath', 'var'), prfOptsPath = [];  end
if ~exist('dataStr', 'var'),     dataStr = 'bold';  end

[session, tasks, runnums] = bidsSpecifyEPIs(projectDir, subject,...
    session, tasks, runnums);


% <dataFolder>
if ~exist('dataFolder', 'var') || isempty(dataFolder)
    dataFolder = 'fmriprep';
end
dataPath = fullfile (projectDir,'derivatives', dataFolder,...
    sprintf('sub-%s',subject), sprintf('ses-%s',session));
if exist(fullfile(dataPath, 'func'), 'dir')
    dataPath = fullfile (dataPath, 'func');
end
rawDataPath = fullfile(projectDir, sprintf('sub-%s', subject), ...
    sprintf('ses-%s', session), 'func');

% <apertureFolder>
if ~exist('apertureFolder', 'var'), apertureFolder = []; end
aperturePath = fullfile(projectDir, 'derivatives', 'stim_apertures', ...
    apertureFolder, sprintf('sub-%s',subject), sprintf('ses-%s',session));
if ~exist(aperturePath, 'dir')
    error('Aperture path not found: %s', aperturePath); 
end
     
%% Create analyzePRF inputs

%****** Required inputs to analyzePRF *******************
% < stimulus>
stimulus = getStimulus(aperturePath, tasks, runnums);

% <data>
data = bidsGetPreprocData(dataPath, dataStr, tasks, runnums);

% <tr>
if ~exist('tr', 'var') || isempty(tr)
    tr = bidsGetJSONval(rawDataPath,tasks, runnums, 'RepetitionTime');
    tr = cell2mat(tr);
    if length(unique(tr)) > 1
        disp(unique(tr))
        error(['More than one TR found:' ...
            'GLMdenoise expects all scans to have the same TR.'])
    else
        tr = unique(tr);
    end
end

%****** Optional inputs to analyzePRF *******************
% prf opts

%% NYI 
% opt = getPRFOpts(prfOptsPath);
opt = [];


%% Run the analyzePRF alogithm
results  = analyzePRF (stimulus,data,tr,opt);

%   <figuredir>
resultsdir   = fullfile (projectDir,'derivatives','analyzePRF', ...
                 sprintf('sub-%s',subject), sprintf('ses-%s',session));

if ~exist(resultsdir, 'dir'); mkdir(resultsdir); end

% save the results
fname = sprintf('sub-%s_ses-%s_results', subject, session);
save(fullfile(resultsdir, fname), 'results', '-v7.3');

end


%% ******************************
% ******** SUBROUTINES **********
% *******************************
function stimulus = getStimulus(aperturePath, tasks, runnums)
% <stimulus> provides the apertures as a cell vector of R x C x time.
%   values should be in [0,1].  the number of time points can differ across runs.

stimFiles = dir(fullfile(aperturePath, '*aperture.mat'));

stimNames = {stimFiles.name};

scan = 1;

for ii = 1:length(tasks)
    for jj = 1:length(runnums{ii})
        
        % check for both 0-padded and non-0-padded integers, throw
        % an error if we find other than 1.
        prfIdx = contains(stimNames,sprintf('task-%s_run-%d_',...
            tasks{ii}, runnums{ii}(jj)));
        prfIdx = or(prfIdx, contains(stimNames,sprintf('task-%s_run-%02d_',...
            tasks{ii}, runnums{ii}(jj))));
        if sum(prfIdx) > 1
           error(['Found more than one stim file for task %s, run %d / '...
                  '%02d'], tasks{ii}, runnums{ii}(jj), runnums{ii}(jj))
        elseif sum(prfIdx) == 0 
            error (['Stim file *_task-%s_run-%d_aperture.mat'...
                    ' was not found'],tasks{ii}, runnums{ii}(jj))
        else
            tmp = load(fullfile(stimFiles(prfIdx).folder,stimFiles(prfIdx).name), 'stimulus');
            stimulus{scan} = tmp.stimulus;
            scan         = scan+1;
        end
    end
end
end


  
           