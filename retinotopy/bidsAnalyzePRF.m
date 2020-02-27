function results = bidsAnalyzePRF(projectDir, subject, session, tasks, runnums, ...
        dataFolder, dataStr, apertureFolder, modelType, prfOptsPath, tr)
%
% results = bidsAnalyzePRF(projectDir, subject, [session], [tasks], [runnums], ...
%        [dataFolder], [apertureFolder], [modelType], [glmOptsPath], [tr]);
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
%     modelType:        name of folder and file to store outputs of analyzePRF (string)
%                           default = apertureFolder;
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
%     modelType         = [];
%
%     % make the stimulus apertures 
%     bidsStimulustoApertures(projectDir, subject, session, tasks, runnums, apertureFolder);
%
%     % run the prf analysis
%     bidsAnalyzePRF(projectDir, subject, session, tasks, runnums, ...
%        dataFolder, dataStr, apertureFolder, modelType, prfOptsPath, tr)
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
 
% <modelType>
if ~exist('modelType', 'var') || isempty(modelType)
    modelType = apertureFolder;
end

%% Create analyzePRF inputs

%****** Required inputs to analyzePRF *******************
% < stimulus>
[stimulus, stimwidthpix] = getStimulus(aperturePath, tasks, runnums);

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

%% Optional inputs to analyzePRF 
%
[averageScans,stimwidthdeg,opt] = getPRFOpts(prfOptsPath);

if ~isempty(averageScans)
   
    % average the requested scans and reduce the stimuli to the unique
    % scans
    [~, ia] = unique(averageScans);
    dims = ndims(data{1});
    
    mn = cell(1,length(ia));
    
    for ii = 1:length(ia)
        whichscans = averageScans==ia(ii);
        datatmp = catcell(dims+1, data(whichscans));
        mn{ii} =  mean(datatmp, dims+1);
    end
    data = mn;
    stimulus = stimulus(ia);
end

%% Save input arguments

inputVar = struct('projectDir', projectDir, 'subject', subject, ...
    'session', session, 'tasks', tasks, 'runnums', runnums, ...
    'dataFolder', dataFolder, 'dataStr', dataStr, 'apertureFolder', apertureFolder, ...
    'modelType', modelType, 'prfOptsPath', prfOptsPath, 'tr', tr, 'stimwidthdeg', stimwidthdeg, ...
    'stimwidthpix', stimwidthpix);
    
fname = sprintf('sub-%s_ses-%s_%s_inputVar.json', subject, session, modelType);

%   <resultsdir>
resultsdir   = fullfile (projectDir,'derivatives','analyzePRF', modelType, ...
                 sprintf('sub-%s',subject), sprintf('ses-%s',session));

if ~exist(resultsdir, 'dir'); mkdir(resultsdir); end

savejson('',inputVar,fullfile(resultsdir,fname));


%% Run the analyzePRF alogithm
results  = analyzePRF (stimulus,data,tr,opt);


% save the results as matlab file
fname = sprintf('sub-%s_ses-%s_%s_results', subject, session, modelType);
save(fullfile(resultsdir, fname), 'results', '-v7.3');

% save the results as mgz files
aPRF2Maps(projectDir, subject, session, modelType);

% save out .png files of angle, ecc, sigma, R2 for lh and rh
Maps2PNG(projectDir, subject, session, modelType);

end


%% ******************************
% ******** SUBROUTINES **********
% *******************************

function [averageScans, stimwidth, opt] = getPRFOpts(prfOptsPath)

if ~exist('prfOptsPath', 'var') || isempty(prfOptsPath)
    prfOptsPath = prfOptsMakeDefaultFile; 
end

json = jsondecode(fileread(prfOptsPath));


if isfield(json, 'averageScans'), averageScans = json.averageScans;
else, averageScans = []; end

if isfield(json, 'stimwidth'), stimwidth = json.stimwidth;
else, error('Stim width not specified in options file'); end

if isfield(json, 'opt'), opt = json.opt; else, opt = []; end

end


function pth = prfOptsMakeDefaultFile()
    % see analyzePRF for descriptions of optional input 
    
    % average scans with identical stimuli
    json.averageScans = [];  % 
    json.stimwidth    = 24.2;  % degrees
    
    % other opts
    json.opt.vxs            = []; 
    json.opt.wantglmdenoise = []; 
    json.opt.hrf            = [];
    json.opt.maxpolydeg     = [];
    json.opt.numperjob      = [];
    json.opt.xvalmode       = [];
    json.opt.seedmode       = [];
    json.opt.maxiter        = [];
    json.opt.display        = 'off';
    json.opt.typicalgain    = [];

                  
    pth = fullfile(tempdir, 'prfOpts.json');
    savejson('', json, 'FileName', pth);
end      
  
           