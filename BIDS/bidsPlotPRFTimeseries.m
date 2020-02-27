function bidsPlotPRFTimeseries(projectDir, subject, session, tasks, runnums, ...
                                dataFolder, dataStr, modelType)

% Example
% projectDir        = '/Volumes/server/Projects/SampleData/BIDS';
% subject           = 'wlsubj042';
% session           = '01';
% tasks             = 'prf';
% runnums           = [1, 2];
% modelType         = 'Fine';
% dataFolder        = 'fmriprep';
% dataStr           = 'fsnative*.mgz';
% 
% bidsPlotPRFTimeseries(projectDir, subject, session, tasks, runnums, ...
%     dataFolder, dataStr, modelType)


%% 0. Define paths and filenames
dataPath             = fullfile(projectDir,'derivatives', dataFolder,...
                        sprintf('sub-%s',subject), sprintf('ses-%s',session), 'func');

resultsFileName      = sprintf('sub-%s_ses-%s_%s_results.mat', ...
                        subject, session, lower(modelType));
                    
aperturePath        = fullfile(projectDir, 'derivatives', 'stim_apertures', sprintf('sub-%s', subject), sprintf('ses-%s', session));   

saveImagesPath      = fullfile(projectDir, 'derivatives', 'analyzePRF', lower(modelType), sprintf('sub-%s', subject), sprintf('ses-%s', session), 'figures');


%% 1. Load results, data, stim and opts

% Load analyzePRF results
load(fullfile(projectDir, 'derivatives', 'analyzePRF', lower(modelType), sprintf('sub-%s', subject), sprintf('ses-%s', session), resultsFileName), 'results');

% Load stim apertures
[stimulus, stimwidthpix] = getStimulus(aperturePath, {tasks}, {runnums});

% Load options use in analyzePRF
prfOpts = loadjson(fullfile(projectDir, 'derivatives', 'stim_apertures', sprintf('sub-%s', subject), sprintf('ses-%s', session), sprintf('prfOpts%s.json', modelType)));

% Load preprocessed data used in analyzePRF
data = bidsGetPreprocData(dataPath, dataStr, {tasks}, {runnums});

% Define variables
res   = [stimwidthpix stimwidthpix]; % row x column resolution of the stimuli
hrf   = results.options.hrf;         % HRF that was used in the model



%% Check if we averaged scans before running analyzePRF
if numel(unique(prfOpts.averageScans)) < length(prfOpts.averageScans)    
    [~,idx] = find(prfOpts.averageScans);   
    averageData = squeeze(nanmean(cat(1,data{idx}),1));
    data = {averageData};
end

%% Prepare the stimuli for use in the model
if iscell(stimulus) && ~exist('averageData', 'var')
    
    % if different stimulus sequences were used (i.e. rings, wedges, bar) then do analysis per stim type
    stimulusPP = {};
    numStimSeq = length(stimulus);
    
    for p = 1:numStimSeq
        stimulusPP{p} = squish(stimulus{p},2)';  % this flattens the image so that the dimensionality is now frames x pixels
        stimulusPP{p} = [stimulusPP{p} p*ones(size(stimulusPP{p},1),1)];  % this adds a dummy column to indicate run breaks
    end
    
else
    stimulusPP = {};
    stimulusPP{1} = squish(stimulus{unique(prfOpts.averageScans)},2)';  % this flattens the image so that the dimensionality is now frames x pixels
    stimulusPP{1} = [stimulusPP{1} ones(size(stimulusPP{1},1),1)];  % this adds a dummy column to indicate run breaks
end


%% Which vertex should we inspect? 

% NOT READY YET: load in ROIs from Benson and Wang Atlases
% [~, ~, bh, bensonAreaLabels] = roisFromAtlas(subject, projectDir, isfmriprep);

[val, vertex] = sort(results.R2, 'descend');
nanIdx = isnan(val);
vertex = vertex(~nanIdx); val = val(~nanIdx);

verticesToPlot = vertex(1:10);

%% For each run, collect the data and the model fit and plot it
datats = {};
modelts = {};
for p= 1:length(data)
    for vx = vertex(1:10)
        
        [x,y] = pol2cart(results.ang(vx),results.ecc(vx));
        prfsz = results.rfsize(vx);
        gain  = results.gain(vx);
        expt  = results.expt(vx);
        
        datats  = data{p}(vx,:)';
        modelts = prf2DGauss(x,y,prfsz,gain,expt,stimulusPP{p}, res, hrf);
        
        % Visualize the results
        figure; clf; hold on;
        set(gcf,'Units','points','Position',[100 100 1000 100]);
        plot(cat(1,datats),'r-');
        plot(cat(1,modelts),'b-');
        xlabel('Time (s)');
        ylabel('BOLD signal');
        title(sprintf('Vertex %d - R2: %1.2f', results.R2(vx)));
        ax = axis;
        axis([0 length(datats) ax(3:4)]);
        title(sprintf('Time-series data - vertex %d - R2: %1.2f', vx, results.R2(vx)));
        
        print(fullfile(saveImagesPath,sprintf('sub-%s_vertex-%d_timeseries',...
            subject, vx)),'-dpng')
    end  
end


end


function modelts = prf2DGauss(x,y,prfsz,gain,expt,stimFrames, res, hrf)
% Define the model function.  This function takes parameters and stimuli as input and
% returns a predicted time-series as output.  Specifically, the variable <pp> is a vector
% of parameter values (1 x 5) and the variable <dd> is a matrix with the stimuli (frames x pixels).
% Although it looks complex, what the function does is pretty straightforward: construct a
% 2D Gaussian, crop it to <res>, compute the dot-product between the stimuli and the
% Gaussian, raise the result to an exponent, and then convolve the result with the HRF,
% taking care to not bleed over run boundaries


% maximum resolution (along any dimension)
resmx   = max(res);                    

% Pre-compute cache for faster execution
[~,xx,yy] = makegaussian2d(resmx,2,2,2,2);

% Create pRF
pRF = makegaussian2d(resmx,x,y,abs(prfsz),abs(prfsz),xx,yy,0,0) / (2*pi*abs(prfsz)^2);

% Crop to resolution
pRF_cropped = vflatten(placematrix(zeros(res), pRF));

% Dot product of pRF with stimulus frames, then raise to exponent
pRFResponse = (stimFrames*[pRF_cropped; 0]) .^ posrect(expt);

% Convolve pRF response with HRF
modelts = conv2run(posrect(gain) * pRFResponse, hrf, stimFrames(:,prod(res)+1));

end

