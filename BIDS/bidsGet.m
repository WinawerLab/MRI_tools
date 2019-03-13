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
    error('The input string does not contain the requested field. \n\tInput string: %s\n\tRequested field: %s', str, param);
else
    val = vals{idx};
end


