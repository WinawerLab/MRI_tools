function writeMGZ(bidsDir, subjectname, data, saveFolder, fileName)

    fsDir = fullfile(bidsDir, 'derivatives','freesurfer');
    
    if ~isfolder(saveFolder)
        mkdir(saveFolder)
    end

    hh = {'lh','rh'};

    % remove conflicting file with duplicate name from path
    conflictingPaths = which('-all', 'mergestruct');
    index = find(contains(conflictingPaths, 'cvncode')); % identify the one of interest
    addpath(fileparts(conflictingPaths{index}), '-begin');

    % get dimensions
    hSize = getSurfsize(fsDir, subjectname); % get hemi size

    hSizeIdx=[1,hSize(1);hSize(1)+1,sum(hSize)];
    % save mgz
    origfile = fullfile(fsDir, subjectname, 'mri', 'orig.mgz');
    mgz = MRIread(origfile);
    
    for hi=1:numel(hh)
        mgz.vol = data(hSizeIdx(hi,1):hSizeIdx(hi,2), 1);
        savePath = fullfile(saveFolder, sprintf('%s.%s.mgz',hh{hi},fileName));
        MRIwrite(mgz, savePath);
    end

end