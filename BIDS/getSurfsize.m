function f = getSurfsize(fsDir, sourcesubject)

fsdirFROM=sprintf('%s/%s/surf',fsDir,sourcesubject);

surf1file = sprintf('%s/%s.sphere.reg%s',fsdirFROM,'lh');
[surf1.vertices,~] = freesurfer_read_surf_kj(surf1file);

surf2file = sprintf('%s/%s.sphere.reg%s',fsdirFROM,'rh');
[surf2.vertices,~] = freesurfer_read_surf_kj(surf2file);


f(1,1) = size(surf1.vertices,1);
f(2,1) = size(surf2.vertices,1);

end