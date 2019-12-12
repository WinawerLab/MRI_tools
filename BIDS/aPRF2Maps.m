function [] = aPRF2Maps(bidsfolder, subject, session, desc)
% Convert the output of analyze PRF, i.e. results struct, to an MGZ 
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

% Path to analyze PRF results
pth = fullfile(bidsfolder, 'derivatives', 'analyzePRF', desc, ['sub-' subject], ['ses-' session]);

% JSON file with input arguments to analyzePRF
opts = loadjson(fullfile(pth, sprintf('sub-%s_ses-%s_%s_inputVar.json', subject, session, desc)));

% AnalyzePRF results file
load(fullfile(pth, sprintf('sub-%s_ses-%s_%s_results.mat', subject, session, desc)), 'results');

% Freesurfer directory
fspth = fullfile(bidsfolder, 'derivatives', 'freesurfer', ['sub-' subject]);

lcurv = read_curv(fullfile(fspth, 'surf', 'lh.curv'));
rcurv = read_curv(fullfile(fspth, 'surf', 'rh.curv'));
assert(isequal(numel(lcurv) + numel(rcurv), numel(results.ang)), ...
    'The number of vertices in the aprf results and the l&r curv files do not match;');

% what are the files expected by bayesian retinotopy?
% left and right theta (in radians), eccen, size, r2