function Maps2PNG(bidsfolder, subject, session, desc)
% Visualize an MGZ ret data on a freesurfer surface in Matlab
%
% Maps2PNG(bidsfolder, subject, session, desc)
%
% INPUTS
%   bidsfolder     : path to BIDS project
%   subject        : BIDS subject name
%   session        : BIDS session name
%   desc           : name of subfolder reflecting model/design that was
%                    used for analysis;

% Example
% bidsfolder = '/Volumes/server/Projects/SampleData/BIDS/';
% subject    = 'wlsubj042';
% session    = '01';
% desc       = 'coarse';



% set up our paths
pth = fullfile(bidsfolder, 'derivatives/analyzePRF/', desc, ['sub-' subject], ['ses-' session]);
freesurfer_dir =  fullfile(bidsfolder, 'derivatives/freesurfer');
setenv('SUBJECTS_DIR',freesurfer_dir)

fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer',  ['sub-' subject], 'surf');
figureDir = fullfile(pth, 'figures');


hemispheres = {'lh';'rh'};
path2roi = {'V1_exvivo';'V2_exvivo'};

% if the fig dir does not exist, make one
if ~exist(figureDir, 'dir')
    mkdir(figureDir);
end


% find the prf data, assuming we want images of everything
mapsList = {'angle', 'eccen', 'sigma', 'vexpl'};
map_file = [];

% load all the data

for hemi = 1 : length(hemispheres)
    
    for thisMap = 1:length(mapsList)
        
        map_file.(mapsList{thisMap}).(hemispheres{hemi}) = load_mgh(fullfile(pth, sprintf('%s.%s.mgz',hemispheres{hemi},mapsList{thisMap})));
        
    end
    
end

% loop through the maps and create png pics


for thisFig = 1:length(mapsList)
    
    figure(1);clf
    for hemi = 1 : length(hemispheres)
        
        roi = [];
        
        for r = 1 : length(path2roi)
            
            ind  = read_label(['sub-' subject],sprintf ('%s.%s%s',hemispheres{hemi},path2roi{r}));
            roi  = [roi; ind(:,1) + 1];
            
        end
        
        curv = read_curv(sprintf('%s/%s.curv',fspth,hemispheres{hemi}));
        surf_file = fullfile(fspth, sprintf('%s.inflated',hemispheres{hemi}));
        [vertices, faces] = read_surf(surf_file);
        
        myroi = zeros(size(map_file.vexpl.(hemispheres{hemi})));
        myroi(roi) = 1;
        
        if contains(desc,'coarse') && thisFig == 1
            
            if nanmean(map_file.vexpl.(hemispheres{hemi})) < 0.11
                map_file.vexpl.(hemispheres{hemi}) = double(map_file.vexpl.(hemispheres{hemi})*100);
            end
            
            thr  = double(map_file.vexpl.(hemispheres{hemi})>0.15) & double(map_file.eccen.(hemispheres{hemi})<10) & myroi ;

        else
            thr  = double(map_file.vexpl.(hemispheres{hemi})>0.15) & double(map_file.eccen.(hemispheres{hemi})<10) & myroi ;
        end
        
        subplot(2,2,hemi)
        plot_mesh(faces, vertices, curv, 'gray');
        freezeColors
        mymap = map_file.(mapsList{thisFig}).(hemispheres{hemi}).*thr;
        mymap(mymap==0) = NaN;
        
        variable = mapsList{thisFig};
        
        switch variable
            
            case 'angle'
                
                plot_mesh(faces, vertices, mymap,'hsv');
                caxis([0 2*pi])
                
            case 'eccen'
                
                plot_mesh(faces, vertices, mymap, 'jet');
                caxis([0 12])
                
            case 'sigma'
                
                plot_mesh(faces, vertices, mymap, 'parula');
                caxis([0 3])
                
                
            case 'vexpl'
                
                plot_mesh(faces, vertices, mymap, 'hot');
                caxis([0 1])
                
        end
        
        if hemi == 1
            set_view(gcf)
            view(45,0)
            
        else
            set_view(gcf)
            view(-45,0)
            
        end
        
    end
    
    subplot(2,2,3:4)
    
    switch variable
        
        case 'angle'
            
            cbh=colorbar('North');
            colormap(hsv)
            caxis([0 2*pi])
            title(variable)
            set(cbh,'XTick',[0 pi 2*pi])
            set(cbh,'XTickLabel',{'0';'\pi';'2\pi'})
            
        case 'eccen'
            
            cbh=colorbar('North');
            colormap(jet)
            caxis([0 12])
            title(variable)
            set(cbh,'XTick',[0 4 8 12])
            set(cbh,'XTickLabel',{'0';'4';'8';'12'})
            
        case 'sigma'
            
            cbh=colorbar('North');
            colormap(parula)
            caxis([0 3])
            title(variable)
            set(cbh,'XTick',[0 1 2 3])
            set(cbh,'XTickLabel',{'0';'1';'2';'3'})
            
        case 'vexpl'
            
            cbh=colorbar('North');
            colormap(hot)
            caxis([0 1])
            set(cbh,'XTick',[0 0.25 0.5 0.75 1])
            set(cbh,'XTickLabel',{'0';'0.25';'0.5';'0.75';'1'})
            
            if contains(desc,'coarse')
                
                title('R2')
            else
                title(variable)
                
            end
    end
    
    
    
    axis off
    set(gca,'Fontsize',30)
    saveas(gcf, ([figureDir (sprintf('/%s_maps.png',  mapsList{thisFig}))]));
    
    
end
%% relevant sub functions
    function plot_mesh(faces, vertices, map, cmap)
        t = trimesh(faces+1, vertices(:,1), vertices(:,2), vertices(:,3), map, 'FaceColor', 'flat');
        t.LineStyle = 'none';
        axis equal; hold on;
        %     cmap = [[1 1 1]*.7; cmap];
        colormap(gcf,cmap);
    end

    function set_view(gcf)
        axis off; set(gcf, 'color','white','InvertHardCopy', 'off');
        view(-45,0);
        material dull;
        h=light; lightangle(h,  45, 45); lighting gouraud;
        h=light; lightangle(h, -45, 45); lighting gouraud;
        h=light; lightangle(h, -45, -90); lighting gouraud;
        set(gcf, 'Position', [150 100 750 625]);
        axis tight
    end

end
