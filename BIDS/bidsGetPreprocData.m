function [data, info] = bidsGetPreprocData(dataPath, tasks, runnums, usePreproc)
%
% Inputs
%   dataPath:   path to folder containing preprocessed data
%   tasks:      BIDS tasks, in cell array
%   runnums:    cell array of runnumbers, equal in length to tasks
%   usePreproc: boolean: if true, use preprocessed data from derivatives
%                       folder, else use data from func folder
%                       default: true
%
% Output
%   data:       the time-series data for each run with dimensions
%                X x Y x Z x time
%   info:       nifti header for each run

if ~exist('usePreproc', 'var') || isempty(usePreproc)
    usePreproc = true;
end

numruns = sum(cellfun(@numel, runnums));

data = cell(1, numruns);
info = cell(1, numruns);

fprintf('Loading data')
scan = 1;
for ii = 1:length(tasks)
    for jj = 1:length(runnums{ii})
        fprintf('.')
        
        if usePreproc
            fnamePrefix  = sprintf('*_task-%s*run-*%d_preproc*',...
                tasks{ii},runnums{ii}(jj));
        else
            fnamePrefix  = sprintf('*_task-%s*run-*%d_bold*',...
                tasks{ii},runnums{ii}(jj));
        end
        
        fname         = dir(fullfile(dataPath, fnamePrefix));
        assert(~isempty(fname));
        
        [~, ~, ext] = fileparts(fname(1).name);
        switch ext
            case {'.nii' '.gz'}
                data{scan}    = niftiread(fullfile (dataPath, fname.name));
                info{scan}    = niftiinfo(fullfile (dataPath, fname.name));
            case '.mgz'
                hemis = {'lh', 'rh'};
                for ll = 1: length (hemis)
                    % We index to make sure the order is always the same
                    idx = contains ({fname.name} , hemis{ll});
                    tempData(ll) = MRIread(fullfile (dataPath, fname(idx).name));
                end
                data{scan}    = cat(2,tempData(1).vol, tempData(2).vol);
                info{scan}    = rmfield(tempData(1), 'vol');
                
            otherwise
                error('Unrecognized file format %s', fname)
        end
        scan          = scan+1;
    end
end
fprintf('\n')
end