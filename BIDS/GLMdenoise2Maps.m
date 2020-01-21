function [] = GLMdenoise2Maps(bidsfolder, subject, session, desc, threshR2beta)
% Convert the output of GLMdenoise, i.e. results struct, to an MGZ, written
% into the GLMdenoise output folder
%
% GLMdenoise2Maps(bidsfolder, subject, session, desc, [threshR2])
% 
% INPUTS
%   bidsfolder     : path to BIDS project
%   subject        : BIDS subject name
%   session        : BIDS session name
%   desc           : type of model [default = ''];
%   threshR2beta   : (optional) explained variance (%) threshold for beta
%                   coefficient maps [default = 2].
%
% OUTPUTS
%
% Example
%  bidsfolder = '/Volumes/server/Projects/BAIR/Data/BIDS/motor/'
%  subject    = 'som756'
%  session    = 'nyu3t01';
%  desc       = 'boldhandandsat';
%  GLMdenoise2Maps(bidsfolder, subject, session, desc)
%
% QQQ why do we get an mv ownership error message?
% Iris Groen Jan 2020

if ~exist('threshR2beta', 'var') || isempty(threshR2beta)
    threshR2beta = 2;
end

% Path to GLMdenoise results
pth = fullfile(bidsfolder, 'derivatives', 'GLMdenoise', desc, ['sub-' subject], ['ses-' session], 'figures');

% GLMdenoise results file
% load(fullfile(pth, sprintf('sub-%s_ses-%s_%s_results.mat', subject, session, desc)), 'results');
d = dir(fullfile(pth, '*results.mat'));
if length(d)>1, warning('found multiple GLMdenoise results files. Loading first one'); end
load(fullfile(d(1).folder,d(1).name),  'results');

% Freesurfer directory
fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer', ['sub-' subject]);

lcurv = read_curv(fullfile(fspth, 'surf', 'lh.curv'));
rcurv = read_curv(fullfile(fspth, 'surf', 'rh.curv'));
assert(isequal(numel(lcurv) + numel(rcurv), numel(results.R2)), ...
    'The number of vertices in the aprf results and the l&r curv files do not match;');

mgz = MRIread(fullfile(fspth, 'mri', 'orig.mgz'));

% left and right hemisphere indices
leftidx  = 1:numel(lcurv);
rightidx = (1:numel(rcurv))+numel(lcurv);

% write out an .mgz with R2 values for each hemisphere
mgz.vol = results.R2(leftidx);
MRIwrite(mgz, fullfile(pth, 'lh.R2.mgz'));
mgz.vol = results.R2(rightidx);
MRIwrite(mgz, fullfile(pth, 'rh.R2.mgz'));

% write out an .mgz with R2 values for each hemisphere and run
n_run = size(results.inputs.design,2);

for ii = 1:n_run    
    % left hemisphere
    mgz.vol = results.R2run(:,leftidx,:,ii);
    MRIwrite(mgz, fullfile(pth, sprintf('lh.R2run%d.mgz', ii)));
    % right hemisphere
    mgz.vol = results.R2run(:,rightidx,:,ii);
    MRIwrite(mgz, fullfile(pth, sprintf('rh.R2run%d.mgz', ii)));
end

% write out an .mgz with the average beta across conditions
mgz.vol = mean(results.modelmd{2}(:,leftidx,:,:),4);
MRIwrite(mgz, fullfile(pth, 'lh.stimcond_mean.mgz'));
mgz.vol = mean(results.modelmd{2}(:,rightidx,:,:),4);
MRIwrite(mgz, fullfile(pth, 'rh.stimcond_mean.mgz'));

% write out an .mgz with the average beta across conditions, thresholded on
% R2
thresh_inx = results.R2<threshR2beta;
tmp = results.modelmd{2};
tmp(:,thresh_inx,:,:) = 0;

mgz.vol = mean(tmp(:,leftidx,:,:),4);
MRIwrite(mgz, fullfile(pth, 'lh.stimcond_mean_threshonR2.mgz'));
mgz.vol = mean(tmp(:,rightidx,:,:),4);
MRIwrite(mgz, fullfile(pth, 'rh.stimcond_mean_threshonR2.mgz'));

% write out an .mgz with betas for each hemisphere and condition
n_cond = size(results.inputs.design{1},2);

for ii = 1:n_cond
    % left hemisphere
    mgz.vol = results.modelmd{2}(:,leftidx,:,ii);
    MRIwrite(mgz, fullfile(pth, sprintf('lh.stimcond%d.mgz', ii)));
    % right hemisphere
    mgz.vol = results.modelmd{2}(:,rightidx,:,ii);
    MRIwrite(mgz, fullfile(pth, sprintf('rh.stimcond%d.mgz', ii)));
end


end


