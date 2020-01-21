function bidsPlotMap(bidsfolder, subject, session, analysis, desc, mapname, hemis, surftype, cmap)
% Visualize an MGZ map as a colorsmale on a freesurfer surface in Matlab 
%
% bidsPlotMap(bidsfolder, subject, session, analysis, desc, mapname, [hemis], [surftype], [cmap])
% 
% INPUTS
%   bidsfolder     : path to BIDS project
%   subject        : BIDS subject name
%   session        : BIDS session name
%   analysis       : BIDS derivatives folder name (e.g. GLMdenoise)
%   desc           : name of subfolder reflecting model/design that was
%                    used for analysis;
%   mapname        : name of the mgz file to dis
%   hemis          : (optional) which hemispheres to plot: options are
%                    'both', 'left', 'right' [default = 'both']
%   surftype       : (optional) freesurfer surface type to be plotted, e.g.
%                    'inflated' [default = 'pial'];
%   cmap           : (optional) colormap used to map values in MGZ file to
%                     a colorscale. Should be an 3-column wide RGB array
%                     [default = 64 length built in Matlab jet colormap].
%
% Note: the colormap is directly mapped to span the range of values in the
% MGZ, to get a meaningful figure it is probably necessary to scale the
% colorscale to a range of interest, for example using the caxis command
% (see example below).
%
% OUTPUTS
%
% Example
%  bidsfolder = '/Volumes/server/Projects/BAIR/Data/BIDS/motor/'
%  subject    = 'som756'
%  session    = 'nyu3t01';
%  analysis   = 'GLMdenoise';
%  desc       = 'boldhandandsat';
%  mapname    = 'stimcond_mean_threshonR2';
%  hemis      = 'right';
%   bidsPlotMap(projectDir, subject, session, analysis, desc, mapname, hemis, 'pial');
%   caxis([0 4]);

% Iris Groen Jan 2020

if ~exist('hemis','var') || isempty(hemis)
  hemis = 'both';
end

if ~exist('surftype','var') || isempty(surftype)
  surftype = 'pial';
end

if ~exist('cmap','var') || isempty(cmap)
  cmap = jet(64);
end

pth = fullfile(bidsfolder, 'derivatives', analysis, desc, ['sub-' subject], ['ses-' session], 'figures');
fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer',  ['sub-' subject], 'surf');

% load vertices and faces 

surf_file_rh = fullfile(fspth, sprintf('rh.%s', surftype));
surf_file_lh = fullfile(fspth, sprintf('lh.%s', surftype));

if exist(surf_file_rh, 'file') && exist(surf_file_lh, 'file')
    fprintf('[%s] Reading Freesurfer surface file %s ...\n',mfilename, surf_file_rh);
    [vertices_r, faces_r] = read_surf(surf_file_rh);
    fprintf('[%s] Reading Freesurfer surface file %s ...\n',mfilename, surf_file_lh);
    [vertices_l, faces_l] = read_surf(surf_file_lh);
else
    fprintf('[%s] Freesurfer surfaces not found - exiting \n',mfilename);
    return
end

% load maps
map_file_rh = fullfile(pth, sprintf('rh.%s.mgz', mapname));
map_file_lh = fullfile(pth, sprintf('lh.%s.mgz', mapname));

if exist(map_file_rh, 'file') && exist(map_file_lh, 'file')
    [map_rh] = load_mgh(map_file_rh);
    [map_lh] = load_mgh(map_file_lh);
else
    fprintf('[%s] No map files found for %s \n',mfilename, mapname);
    return
end

%set nans to zero otherwise those vertices won/t be visible
map_rh(isnan(map_rh)) = 0;
map_lh(isnan(map_lh)) = 0;

% Plot mesh
figure;hold on

switch hemis
    case 'both'
        plot_mesh(faces_r, vertices_r, map_rh, cmap);
        plot_mesh(faces_l, vertices_l, map_lh, cmap);
	case 'right'
        plot_mesh(faces_r, vertices_r, map_rh, cmap);
    case 'left'
        plot_mesh(faces_l, vertices_l, map_lh, cmap);
end

% Set view parameters
set_view(gcf)

% Add colorbar
colorbar;

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