function val = bidsGetJSONval(rawDataPath,tasks, runnums, fieldname)
% Return the value(s) of a field in json files associated with functional
% data in a BIDS directory
%
% Inputs:
%   rawDataPath: 
%   tasks:
%   runnums: 
%   fieldname:  
%
%  Output
%   val:
%
% Example: 

val = cell(1,length(tasks));

scan = 1;
for ii = 1:length(tasks)
    for jj = 1:length(runnums{ii})
        
        jsonPrefix        = sprintf('*_task-%s*run-%d_bold.json',tasks{ii}, runnums{ii}(jj));
        jsonPrefixZeroPad = sprintf('*_task-%s*run-%02d_bold.json',tasks{ii}, runnums{ii}(jj));
        jsonName          = dir(fullfile(rawDataPath, jsonPrefix));
        if  ~strcmp(jsonName,jsonPrefixZeroPad)
            jsonName    = [jsonName ;dir(fullfile(rawDataPath, jsonPrefixZeroPad))];
        end
         % This guarantees that we found at least one
        assert(~isempty(fname));
        
        json        = fileread(fullfile (rawDataPath, jsonName.name));
        jsonInfo    = jsondecode(json);
        val{scan}   = jsonInfo.(fieldname);
        scan        = scan+1;
    end
end

end
