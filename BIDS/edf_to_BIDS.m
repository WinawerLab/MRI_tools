%% edf_to_BIDS

% Convert EDF eyetracking files to BIDS compatible _physio files
% For example: 
%  sub-control01_task-nback_physio.tsv.gz
%  sub-control01_task-nback_physio.json

% Dependencies:
% mgl and mrTools toolboxes
% addpath(genpath('/Volumes/server/Projects/MEG/Eyetracking_scripts/toolboxes'))

% jsonwrite
% tbUse('JSONio')


%% 1. Define parameters and paths

% BIDS specifics
sub  = 'wlsubj042';
ses  = 'mri3t01';
task = 'sfp';

% Get edf path
rawDataPath = '/Volumes/server/Projects/tutorial_dataset/spatial_frequency_preferences/raw_behavior';
savePath    = rawDataPath;

%% 1. Find edf data

% Find EDF file(s)
edffile = dir(fullfile(rawDataPath,'*.edf'));
numRuns = size(edffile,1);

%% 1. Load edf data

for r = 1:numRuns
    % EYD is a struct that contains the data and additional info in the edf
    % file
    eyd = mglEyelinkEDFRead(fullfile(rawDataPath,edffile(r).name));

    % Convert struct to cell
    dataCells = struct2cell(eyd.gaze);


    % Define BIDS compatible filename
    tsvFileName = sprintf('sub-%s_ses-%s_task-%s_run-0%d_physio.tsv.gz', sub,ses,task,r);

    % Write data cell to a tsv file.
    dlmwrite(fullfile(savePath,tsvFileName),dataCells, 'Delimiter', '\t');

end

return
%% Create a JSON file

jsonFileName = sprintf('sub-%s_ses-%s_task-%s_physio.json', sub,ses,task);


% Question: do we want to add more meta Data?
S = struct( ...
        'SamplingFrequency', eyd.samplerate, ... % Hz
        'FileFormat', 'EDF', ...
        'StartTime', eyd.gaze.time(1), ...
        'Columns', fieldnames(eyd.gaze), ...
        'RecordedEye', eyd.whicheye);
    
% Save physio json file    
jsonwrite(fullfile(savePath, jsonFileName),S);

