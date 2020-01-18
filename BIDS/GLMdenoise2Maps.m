function [] = GLMdenoise2Maps(bidsfolder, subject, session, desc)
% Convert the output of GLMdenoise, i.e. results struct, to an MGZ 
%
%
% INPUTS
%   bidsfolder  : path to BIDS project
%   subject     : BIDS subject name
%   session     : BIDS session name
%   desc        : type of model [default = ''];
%
% OUTPUTS
%
% Example
%  bidsfolder =  '/Volumes/server/Projects/SampleData/BIDS/'
%  subject    =  'wlsubj042'
%  session    = '01';
%  desc       = 'coarse';
%  GLMdenoise2Maps(bidsfolder, subject, session, desc)

% Path to GLMdenoise results
pth = fullfile(bidsfolder, 'derivatives', 'GLMdenoise', desc, ['sub-' subject], ['ses-' session], 'figures');

% GLMdenoise results file
load(fullfile(pth, sprintf('sub-%s_ses-%s_%s_results.mat', subject, session, desc)), 'results');

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
    mgz.vol = squeeze(results.R2run(:,leftidx,:,ii));
    MRIwrite(mgz, fullfile(pth, sprintf('lh.R2run%d.mgz', ii)));
    % right hemisphere
    mgz.vol = squeeze(results.R2run(:,rightidx,:,ii));
    MRIwrite(mgz, fullfile(pth, sprintf('rh.R2run%d.mgz', ii)));
end

% write out an .mgz with betas for each hemisphere and condition
n_cond = size(results.inputs.design{1},2);
for ii = 1:n_cond
    % left hemisphere
    mgz.vol = squeeze(results.modelmd{2}(:,leftidx,:,ii));
    MRIwrite(mgz, fullfile(pth, sprintf('lh.stimcond%d.mgz', ii)));
    % right hemisphere
    mgz.vol = squeeze(results.modelmd{2}(:,rightidx,:,ii));
    MRIwrite(mgz, fullfile(pth, sprintf('rh.stimcond%d.mgz', ii)));
end

% QQQ what to do about the mv ownership error message?

end


