%% adjust this to your situation>
addpath('~/scripts/neuralAdapatationInMTL'); % scripts in this repo: https://github.com/rebrowski/neuralAdapatationInMTL
stimdir = '~/scripts/neuralAdapatationInMTL/stimuli/'; % folder with images contained in this repo

% path to large datafiles that need to be downloaded separately 
% like the iEEG traces (ieeg_linked_mastoids_256Hz.mat) and single unit data (sessions.mat)
if isunix
    datadir = '/media/treber/neuralAdaptationInMTL';
end

% isunix returns "true" also for mac, so...
if ismac
    datadir = '~/OneDrive/data/neuralAdaptationInMTL/'; 
end

P_W_D = pwd;
cd('~/scripts/IOSRToolbox');% install once from: https://github.com/IoSR-Surrey/MatlabToolbox
+iosr.install;
cd(P_W_D);
