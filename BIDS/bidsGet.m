function val = bidsGet(str, param, varargin)
%
% Example:
%   str = 'sub-wlsubj050_ses-nyu3T01_task-temporalpattern_run-3_preproc.nii.gz';
%   param = 'ses';
%   val = bidsGet(str, param);
%
% what is the str is a path rather than a filename?
%   for example
%   str = '/Volumes/server/Projects/BAIR/Analyses/visual/sub-wlsubj050/ses-nyu3t01'
%   param = 'ses';
%   val = bidsGet(str, param);


splitstrings = strsplit(str, {'_' filesep '.'},  'CollapseDelimiters',true);
    

for ii = 1:length(splitstrings)
   tmp = regexp(splitstrings{ii}, '-', 'split');
   
   if length(tmp) == 2
        
       params{ii} = tmp{1};
       vals{ii}   = tmp{2};
   end
end


idx = find(contains(params, param));
    
if isempty(idx)
    warning('The input string does not contain the BIDS formatted requested field. \n\tInput string: %s\n\tRequested field: %s', str, param);
    
    % work around to accomodate non-BIDS compliant format task-<taskname><runnumber>
    taskIdx = find(strcmp(params, 'task'));
    if ~isempty(taskIdx)  % Ensure 'task' exists in params
        taskValue = vals{taskIdx};  % Extract the task name
        
        taskName = regexprep(taskValue, '\d+$', '');  % Remove trailing digits
    
        % Display results
        disp(['Original task value: ', taskValue]);
        disp(['Task name without numbers: ', taskName]);
    else
        error('Task parameter not found in params array.');
    end

    val = nonbidsGetrun(str, taskName);

else
    val = vals{idx};
end


