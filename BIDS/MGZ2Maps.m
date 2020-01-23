function MGZ2Maps(bidsfolder, subject, session, desc, cmap)

% Visualize an MGZ ret data on a freesurfer surface in Matlab 
%
% MGZ2Maps(bidsfolder, subject, session, desc, cmap)
% 
% INPUTS
%   bidsfolder     : path to BIDS project
%   subject        : BIDS subject name
%   session        : BIDS session name
%   desc           : name of subfolder reflecting model/design that was
%                    used for analysis;
%   cmap           : (optional) colormap used to map values in MGZ file to
%                     a colorscale. Should be an 3-column wide RGB array
%                     [default = 64 length built in Matlab jet colormap].
%
% Note: the colormap is directly mapped to span the range of values in the
% MGZ, to get a meaningful figure it is probably necessary to scale the
% colorscale to a range of interest, for example using the caxis command
% (see example below).
%
% Example
% bidsfolder = '/Volumes/server/Projects/SampleData/BIDS/';
% subject    = 'wlsubj042';
% session    = '01';
% desc       = 'coarse';

%caxis([0 4]);


if ~exist('cmap','var') || isempty(cmap)
  cmap = jet(64);
end

% set up our paths
pth = fullfile(bidsfolder, 'derivatives/analyzePRF/', desc, ['sub-' subject], ['ses-' session]);
fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer',  ['sub-' subject], 'surf');
figureDir = fullfile(pth, 'figures');

% if the fig dir does not exist, make one
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end

% load vertices and faces 
surf_file_rh = fullfile(fspth, sprintf('rh.inflated'));
surf_file_lh = fullfile(fspth, sprintf('lh.inflated'));

if exist(surf_file_rh, 'file') && exist(surf_file_lh, 'file')
    fprintf('[%s] Reading Freesurfer surface file %s ...\n',mfilename, surf_file_rh);
    [vertices_r, faces_r] = read_surf(surf_file_rh);
    fprintf('[%s] Reading Freesurfer surface file %s ...\n',mfilename, surf_file_lh);
    [vertices_l, faces_l] = read_surf(surf_file_lh);
else
    fprintf('[%s] Freesurfer surfaces not found - exiting \n',mfilename);
    return
end

% find the prf data, assuming we want images of everything
mapsList = {'angle', 'eccen', 'sigma', 'vexpl'};
for thisMap = 1:length(mapsList)
    map_file_lh{thisMap,:} = fullfile(pth, sprintf('lh.%s.mgz', mapsList{1,thisMap}));
    map_file_rh{thisMap,:} = fullfile(pth, sprintf('rh.%s.mgz', mapsList{1,thisMap}));
end

% and load in the maps
for thisMap = 1:length(map_file_lh)
    map_lh{thisMap,:} = load_mgh(map_file_lh{thisMap,1});
    map_rh{thisMap,:} = load_mgh(map_file_rh{thisMap,1});
end

%loop though and set all nans to zero otherwise those vertices will not be visible in the mesh
for thisMap = 1:length(map_lh)
    map_lh{thisMap,:}(isnan(map_lh{thisMap,:})) = 0;
    map_rh{thisMap,:}(isnan(map_rh{thisMap,:})) = 0;
end


% Plot and save out left hemi figs
for thisFig = 1:length(map_lh)
    clf
    plot_mesh(faces_l, vertices_l, map_lh{thisFig,1}, cmap);
    set(gcf,'Visible','off');
    title((sprintf('LH.%s', mapsList{1,thisFig})));
    set_view(gcf)
    saveas(gcf, ([figureDir (sprintf('/LH_%s.png', mapsList{1,thisFig}))]));
end

%and the right hemi
for thisFig = 1:length(map_rh)
    clf
    plot_mesh(faces_r, vertices_r, map_rh{thisFig,1}, cmap);
    set(gcf,'Visible','off');
    title((sprintf('RH.%s', mapsList{1,thisFig})));
    set_view(gcf)
    saveas(gcf, ([figureDir (sprintf('/RH_%s.png', mapsList{1,thisFig}))]));
end

%% relevant sub functions
function plot_mesh(faces, vertices, map, cmap)
    t = trimesh(faces+1, vertices(:,1), vertices(:,2), vertices(:,3), map, 'FaceColor', 'flat'); 
    t.LineStyle = 'none';    
	axis equal; hold on;
    cmap = [[1 1 1]*.7; cmap];
    colormap(gcf,cmap);    
end

function set_view(gcf)
    axis off; set(gcf, 'color','white','InvertHardCopy', 'off');
    view(-45,0);
    material dull;
    h=light; lightangle(h,  45, 45); lighting gouraud;
    h=light; lightangle(h, -45, 45); lighting gouraud;
    h=light; lightangle(h, -45, -90); lighting gouraud;
    set(gcf, 'Position', [150 100 1500 1250]);
    axis tight    
end
end

