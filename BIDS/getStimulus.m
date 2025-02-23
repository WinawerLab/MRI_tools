function [stimulus, stimwidthpix] = getStimulus(aperturePath, tasks, runnums)
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
            tmp = load(fullfile(stimFiles(prfIdx).folder,stimFiles(prfIdx).name)); %, 'stimulus');
            try
                stimulus{scan} = tmp.stimulus;
            catch
                stimulus{scan} = tmp.stim;
            end
            scan         = scan+1;
        end
    end
    
    stimwidthpix = size(stimulus{1},2);

end
end