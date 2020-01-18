function bidsPlotMap(bidsfolder, subject, session, analysis, desc, mapname, surftype, cmap)

% bidsPlotMap(bidsfolder, subject, session, analysis, desc, mapname, [surftype], [cmap])

if ~exist('surftype','var') || isempty(surftype)
  surftype = 'pial';
end

if ~exist('cmap','var') || isempty(cmap)
  cmap = jet(64);
end

pth = fullfile(bidsfolder, 'derivatives', analysis, desc, ['sub-' subject], ['ses-' session], 'figures');
fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer',  ['sub-' subject], 'surf');

surf_file_rh = fullfile(fspth, sprintf('rh.%s', surftype));
surf_file_lh = fullfile(fspth, sprintf('lh.%s', surftype));

map_file_rh = fullfile(pth, sprintf('rh.%s.mgz', mapname));
map_file_lh = fullfile(pth, sprintf('lh.%s.mgz', mapname));

% load files
if exist(surf_file_rh, 'file') && exist(surf_file_lh, 'file')
    fprintf('[%s] Reading Freesurfer surface file %s ...\n',mfilename, surf_file_rh);
    [vertices_r, faces_r] = read_surf(surf_file_rh);
    fprintf('[%s] Reading Freesurfer surface file %s ...\n',mfilename, surf_file_lh);
    [vertices_l, faces_l] = read_surf(surf_file_lh);
else
    fprintf('[%s] Freesurfer surfaces not found - exiting \n',mfilename);
    return
end

if exist(map_file_rh, 'file') && exist(map_file_lh, 'file')
    [map_rh] = load_mgh(map_file_rh);
    [map_lh] = load_mgh(map_file_lh);
else
    fprintf('[%s] No map files found for %s \n',mfilename, mapname);
    return
end

% Plot mesh
figure;hold on
plot_mesh(faces_r, vertices_r, map_rh, cmap);
plot_mesh(faces_l, vertices_l, map_lh, cmap);

% Set view parameters
set_view(gcf)

end

function plot_mesh(faces, vertices, map, cmap)
    t = trimesh(faces+1, vertices(:,1), vertices(:,2), vertices(:,3), map, 'FaceColor', 'flat'); 
    t.LineStyle = 'none';    
	axis equal; hold on;
    cmap = [[1 1 1]*.7; cmap];
    colormap(gcf,cmap);    
end

function set_view(gcf)
    axis off; set(gcf, 'color','white','InvertHardCopy', 'off');
    view(0,0);
    material dull;
    h=light; lightangle(h,  45, 45); lighting gouraud;
    h=light; lightangle(h, -45, 45); lighting gouraud;
    h=light; lightangle(h, -45, -90); lighting gouraud;
    set(gcf, 'Position', [150 100 1500 1250]);
    axis tight    
end