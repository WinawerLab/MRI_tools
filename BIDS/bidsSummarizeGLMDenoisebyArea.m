function [meanBeta,betasSE , GLMconditions] = bidsSummarizeGLMDenoisebyArea (projectDir,isfmriprep,...
    modelType, subject,session,tasks,conditionsOfInterest, makeFigures, saveFigures)
% Computes and plots a mean beta weight for each stimulus condition and
%   visual area using GLM denoise results from surface time-series
%
% summarizeGLMDenoisebyArea (projectDir , subject, modelType, [session],...
%   [tasks], [conditionsOfInterest], [makeFigures], [saveFigures])
%
% Required input:
%
%   projectDir : path where the BIDS projects lies (string)
%   subject    : BIDS subject name (string, all lower case)
%   modelType  : name of folder containing outputs of GLMdenoised (string)
%   isfmriprep : whether project is an fMRIprep project
%
% Optional input:
%
%   session              : BIDS session name (string, all lower case)
%   tasks                : The tasks used for running the GLM
%                               default : use three tasks 
%                                   (temporalpattern, spatialpattern, spatialobject)
%                               Note: The total number of conditions should
%                               match the number of columns in design matrices
%                               used for GLM
%   conditionsOfInterest : one or more types experimental conditions used
%                           for thresholding (string or cell array of strings)
%                               default: uses all conditions
%   makeFigures          : 0 = don't make a plot in this function (default: false)
%   saveFigures          : 0 = don't make a plot in this function (default: false)
%
% example 1:
%
%   projectDir = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%   isfmriprep = true;
%   subject    = 'wlsubj051';
%   modelType  = 'threeTasksUpsampledAssumeHRF';
%
%   bidsSummarizeGLMDenoisebyArea(projectDir,isfmriprep,modelType,subject)
%
% example 2:
%
%   projectDir           = '/Volumes/server/Projects/BAIR/Data/BIDS/visual';
%   isfmriprep           = true;
%   subject              = 'wlsubj062';
%   modelType            = 'twoTasksUpsampledAssumeHRF';
%   conditionsOfInterest = ["crf"];
%   session              = 'nyu3t03';
%   makeFigures          = true;
%   saveFigures          = true;
%   tasks                = {'spatialpattern','spatialobject'};
%
%   bidsSummarizeGLMDenoisebyArea (projectDir,isfmriprep,...
%    modelType, subject,session,tasks,conditionsOfInterest, makeFigures, saveFigures)
%

%% Check inputs and define paths
if ~exist('session', 'var') || isempty(session)
    [session, ~, ~] = bidsSpecifyEPIs(projectDir,subject);
end
if ~exist ('modelType', 'var' )|| isempty(modelType)
    error ('modelType not defined')
end
if ~exist ('isfmriprep', 'var' )|| isempty(isfmriprep)
    error ('isfmriprep not defined')
end
if ~exist ('saveFigures', 'var' )|| isempty(saveFigures)
    saveFigures = 0;
end
if ~exist ('makeFigures', 'var' )|| isempty(makeFigures)
    makeFigures = 0;
end
% set path to GLM results
resultsDir = fullfile(projectDir, 'derivatives','GLMdenoise', modelType,...
    sprintf('sub-%s',subject), sprintf('ses-%s',session), 'figures');

%% Load Data, find stimulus conditions and set ROI labels

% load in ROIs from Benson and Wang Atlases
[~, ~, bh, bensonAreaLabels] = roisFromAtlas(subject, projectDir, isfmriprep);

% load in GLM Denoised results
files = dir(fullfile(resultsDir,'*results.mat'));
load(fullfile(resultsDir, files(1).name), 'results');

% get beta weights and R2 from GLM denoise results
betasmd = squeeze(results.modelmd{2});
betas   = squeeze(results.models{2});
R2      = results.R2';
clear results;

% for now load stimulus condition categories from matfile
load('designMatrixConditions.mat', 'allConditions', 'GLMdefaults','temporalpattern',...
    'spatialobject', 'spatialpattern', 'conditionSubsets');

GLMconditions = [];
if ~exist ('tasks', 'var' )|| isempty(tasks)
    GLMconditions = GLMdefaults;
else
    experiments = {'spatialpattern','temporalpattern', 'spatialobject'};
    %loop through the experiments to make sure the condition order is
    %correct and the same as our design matrices
    for tt = 1:length(experiments)
        if  any(contains(tasks , experiments{tt}, 'IgnoreCase',true))
            tmp = eval(experiments{tt});
            GLMconditions = [GLMconditions tmp];
        else %do nothing
        end
    end
end
 
% Use all conditions, or find conditions that match 'conditionsOfInterest'
if ~exist('conditionsOfInterest', 'var') || isempty(conditionsOfInterest)
    conditionsToUse = ones(1,size(betasmd,2));
else
    conditionsToUse = contains(GLMconditions, conditionsOfInterest, 'IgnoreCase',true);
end

%% Calculate mean and weighted mean beta weights by condition

% Set eccetricity limits (in degrees)
eccenLmt = bh.eccen > 1 & bh.eccen < 8;

% set all negative beta weights to NaN
negativeBetas = mean(betasmd(:,conditionsToUse),2)<0;
betasmd(negativeBetas,:) = nan;
negativeBetasSE = squeeze(mean(betas(:,conditionsToUse,:),2)<0);
betas(negativeBetasSE(1),:,negativeBetasSE(2)) = nan;

% Find all vertices with an R2 less than 2% and set them to NaN
idx = R2 < 2;
betasmd(idx, :) = nan;
betas(idx,:,:)  = nan;

% meanBeta is a matrix of beta weights, n x m, where n is the number of
%   ROIs and m is the number of stimulus conditions

% Loop through visual areas and average beta weights for each condition
for ii = 1:(length(unique(bh.varea))-1)
    indices =  eccenLmt & bh.varea == ii;
    betasSE(ii,:)  = std(squeeze(nanmean(betas(indices, :,:))), [], 2)';
    for jj = 1: length(GLMconditions)
        meanBeta(ii,jj) = nanmean(betasmd(indices, jj));
    end
end

%% Plot it
if makeFigures
    for ll = 1:length(bensonAreaLabels)
        f = figure('Name', bensonAreaLabels{ll});
        
        bar(meanBeta(ll,:)), hold on
        errorbar(meanBeta(ll,:), betasSE(ll,:),'LineStyle','none')
        xticks(1:length(GLMconditions))
        xticklabels([])
        title (bensonAreaLabels{ll})
        ylabel ('BW')
        xticklabels(GLMconditions)
        xtickangle (45)
        set(f, 'Position', (f.Position.*[1 1 1.5 2])) 
    end
    xticklabels(GLMconditions)
    xtickangle (45)
    
    % save it
    if(saveFigures)
        imDir = fullfile(projectDir, 'derivatives','GLMdenoise', modelType);
        print (fullfile(imDir,sprintf('sub-%s_stimSubset-%d_AvgBW',...
            subject{jj}, ii)),'-dpng')
    end
end

% Voxel selection. All analyses were restricted to voxels that satisfy the
% following three criteria. First voxels must be located within 2°?10°
% (eccentricity) based on the pRF model. Second, voxels must have a
% positive Betaweight for the average across all nonblank temporal
% conditions (and averaged across bootstraps). The bootstraps, computed
% byGLMdenoise, were derived by sampling with replacement from the repeated
% scans. Third, voxels must have 2% GLM R2. Voxels that satisfy all
% criteria were averaged within a participant to yield 13 Beta weights per
% ROI per participant per 100 bootstraps. The data were then averaged
% across participants. Averaging within a participant before averaging
% across participants ensures that the contribution for each participant
% has the same weight, regardless of the numbers of voxels per participant
