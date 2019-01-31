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

val = cell(1,length(tasks));

scan = 1;
for ii = 1:length(tasks)
    for jj = 1:length(runnums{ii})
        
        jsonPrefix = sprintf('*_task-%s_run-%d_bold.json',tasks{ii}, runnums{ii}(jj));
        jsonName    = dir(fullfile(rawDataPath, jsonPrefix));
        json        = fileread(fullfile (rawDataPath, jsonName.name));
        jsonInfo    = jsondecode(json);
        val{scan}   = jsonInfo.(fieldname); 
        scan        = scan+1;
    end
end

end