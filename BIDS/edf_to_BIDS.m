%% edf_to_BIDS

% Convert EDF eyetracking files to BIDS compatible _physio files
% For example: 
%  sub-control01_task-nback_physio.tsv.gz
%  sub-control01_task-nback_physio.json

% Dependencies:
% mgl and mrTools toolboxes
% addpath(genpath('/Volumes/server/Projects/MEG/Eyetracking_scripts/toolboxes'))

%% 1. Define parameters and paths

% BIDS specifics
sub  = 'wlsubj042';
ses  = 'mri3t01';
task = 'pRF';

% Get edf path
rawDataPath = '/Volumes/server/Projects/tutorial_dataset/color_pRF/raw_behavior';
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
    bidsFileName = sprintf('sub-%s_ses-%s_task-%s_run-%d_physio.tsv.gz', sub,ses,task,r);

    % Write data cell to a tsv file.
    dlmwrite(fullfile(savePath,bidsFileName),dataCells, 'Delimiter', '\t');

end


%% Create a JSON file


metaData = struct( ...
        'SamplingFrequency', samplerate, ... % Hz
        'FileFormat', 'EDF', ...
        'StartTime', eyd.gaze.time(1), ...
        'Columns', fieldnames(eyd.gaze), ...
        'RecordedEye', eyd.whicheye);
    
    


