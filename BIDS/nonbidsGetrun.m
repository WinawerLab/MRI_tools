function val = nonbidsGetrun(str, task)
%
% Example:
%   str = 'sub-wlsubj050_ses-nyu3T01_task-temporalpattern3_preproc.nii.gz';
%   task = 'temporalpattern';
%   val = nonbidsGetrun(str, task);


splitstrings = strsplit(str, {'_' filesep '.'},  'CollapseDelimiters',true);
    

for ii = 1:length(splitstrings)
   tmp = regexp(splitstrings{ii}, '-', 'split');
   
   if length(tmp) == 2
        
       params{ii} = tmp{1};
       vals{ii}   = tmp{2};
   end
end


taskValue = vals{contains(vals, task)};

% Extract trailing numbers using regexp
idx = regexp(taskValue, '\d+$', 'match');

% Convert from cell array to a string (if found)
if ~isempty(idx)
    idx = idx{1};  % Extract the matched number
else
    idx = '';  % If no numbers found, return an empty string
end
    
if isempty(idx)
    error('The input string does not contain the requested field. \n\tInput string: %s\n\tRequested field: %s', str, task);
else
    val = idx;
end

end


