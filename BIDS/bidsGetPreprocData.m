function [data, info, fullFile] = bidsGetPreprocData(dataPath, dataStr, tasks, runnums)
%
% Inputs
%   dataPath:   path to folder containing preprocessed data
%   dataStr:    text string to specify filename for data 
%   tasks:      BIDS tasks, in cell array
%   runnums:    cell array of runnumbers, equal in length to tasks
%
% Output
%   data:       the time-series data for each run with dimensions
%                X x Y x Z x time
%   info:       nifti header for each run
%   fullFile:   The (nifti) time-series data including header information 
%
% Example: 


numruns = sum(cellfun(@numel, runnums));

data = cell(1, numruns);
info = cell(1, numruns);

fprintf('Loading data')
scan = 1;
for ii = 1:length(tasks)
    for jj = 1:length(runnums{ii})
        fprintf('.')
        
        % we want to check for both 0-padded and non-0-padded versions...
        fnamePrefix  = sprintf('*_task-%s*run-%d_*%s*',...
            tasks{ii},runnums{ii}(jj), dataStr);
        fnamePrefixZeroPad = sprintf('*_task-%s*run-%02d_*%s*',...
            tasks{ii},runnums{ii}(jj), dataStr);
        

        fname = dir(fullfile(dataPath, fnamePrefix));
        % we only need to check both if they're different; if we're
        % looking at run 10, 0-padded and non-0-padded will be the
        % same string
        if ~strcmp(fnamePrefix, fnamePrefixZeroPad)
            fname = [fname; dir(fullfile(dataPath, fnamePrefixZeroPad))];
        end
        
        % We want brain images, not text files, so remove json/tsv files
        istxt = contains({fname.name}, {'.json', '.tsv'});
        fname = fname(~istxt);
             
        % This guarantees that we found at least one
        assert(~isempty(fname));
        
        [~, ~, ext] = fileparts(fname(1).name);
        switch ext
            case {'.nii' '.gz'}
                % if these are niftis, then there should only be
                % one file
                assert(length(fname) == 1);
                fullFile{scan}= niftiRead(fullfile (dataPath, fname.name));
                data{scan}    = fullFile{scan}.data; 
                info{scan}    = niftiinfo(fullfile (dataPath, fname.name));
            case '.mgz'
                % if they're surfaces, there should be two of them
                assert(length(fname) == 2);                
                hemis = {'hemi-L', 'hemi-R'}; % hemis = {'lh', 'rh'};
                for ll = 1: length (hemis)
                    % We index to make sure the order is always the same
                    idx = contains ({fname.name} , hemis{ll}, 'IgnoreCase', true);
                    tempData(ll) = MRIread(fullfile (dataPath, fname(idx).name));
                end
                fullFile{scan}= [];
                data{scan}           = cat(2,tempData(1).vol, tempData(2).vol);
                info{scan}           = rmfield(tempData(1), 'vol');
                info{scan}.ImageSize = size(tempData(1).vol);
            case '.gii'
                assert(length(fname) == 2);
                hemis = {'L', 'R'};
                for ll = 1: length (hemis)
                    % We index to make sure the order is always the same
                    idx = contains ({fname.name} , sprintf('hemi-%s',hemis{ll}));
                    tempData{ll} = gifti(fullfile (dataPath, fname(idx).name));
                end
                fullFile{scan}= [];
                data{scan}    = cat(1,tempData{1}.cdata, tempData{2}.cdata);
            otherwise
                error('Unrecognized file format %s', ext)
        end
        scan          = scan+1;
    end
end
fprintf('\n')