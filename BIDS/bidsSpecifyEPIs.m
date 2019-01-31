function [session, tasks, runnum] = bidsSpecifyEPIs(projectDir, subject,...
    session, tasks, runnum)
% Specify tasks and run numbers and verify paths for BIDS session
% [session, tasks, runnum] = bidsSpecifyEPIs(projectDir, subject, [session], [tasks], [runnum])
%
% Input
%     projectDir:       path where the BIDS projects lies (string)
%     subject:          BIDS subject name (string, all lower case)
%     session:          BIDS session name (string, all lower case)
%                           default: folder name inside subject dir (if
%                           more than one or less than one, return error)
%     tasks:            one or more BIDS tasks (string or cell array of strings)
%                           default: all tasks in session
%     runnum:           BIDS run numbers (vector or cell array of vectors)
%                           default: all runs for specified tasks
% Output
%     session:          String. See input
%     tasks:            Cell array of strings. See input
%     runnum:           Cell array of vectors. See input 
%
% Example 1
%     projectDir        = '/Volumes/server/Projects/BAIR/Data/BIDS/visual'; 
%     subject           = 'wlsubj054';
%     [session, tasks, runnum] = bidsSpecifyEPIs(projectDir, subject)
%
% See also bidsGLM.m bidsTSVtoDesign.m


% <projectDir>
if ~exist('projectDir', 'var') || isempty(projectDir)
    error('projectDir not defined');
end    

% <subject>
if ~exist('subject', 'var') || isempty(subject)
    error('subject not defined');
end


subjectDir = fullfile(projectDir, sprintf('sub-%s', subject));
if ~exist(subjectDir, 'dir')
    error('subject dir not found: %s', subjectDir); 
end

%% Set the optional inputs

% <session>
if ~exist('session', 'var') || isempty(session)
    d = dir(fullfile(subjectDir, 'ses-*'));
    if isempty(d), error('No session folders found in %s', subjectDir); end
    if length(d)>1, error('Two or more session folders found in %s', subjectDir); end
    session = bidsGet(d.name, 'ses');
end
sessionDir = fullfile(subjectDir, sprintf('ses-%s', session));

% <tasks>
if ~exist('tasks', 'var') || isempty(tasks)
    d = dir(fullfile(sessionDir, 'func', '*bold.nii*'));
    
    taskname = cell(1,length(d));
    for ii = 1:length(d)
       taskname{ii} = bidsGet(d(ii).name, 'task');
    end
    tasks = unique(taskname);       
end
if ~iscell(tasks), tasks = {tasks}; end

% <runnum>
if ~exist('runnum', 'var') || isempty(runnum)
         
    runnum = cell(1,length(tasks));
    for ii = 1:length(tasks)
       files = dir(fullfile(sessionDir, 'func', sprintf('*task-%s_*bold.nii.*', tasks{ii})));
       for jj = 1:length(files)
           runnum{ii}(jj) = str2double(bidsGet(files(jj).name, 'run'));
       end       
    end
          
end
if ~iscell(runnum), runnum = {runnum}; end

end