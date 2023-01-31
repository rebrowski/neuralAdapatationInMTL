clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

min_n_re_stim = 4; % inclusion criterion for units responding to at least min_n_re_stimuli
load(sprintf('tuningCurvesMin%dResponsesPerUnit.mat', min_n_re_stim));

%% plotting params
xl = [0.5, min_n_re_stim]; % how many responses per unit are considered
xt = [1:1:xl(2)];
yl = [0.6 1];
yt = [0.5:0.1:1];
alphacrit = 0.05;
    
yl = [0.5 1.1];
figh = figure('Color', 'w');
figh.PaperUnits = 'inches';
figh.PaperPosition = [0 0 7.4, 3.6];

subplot(1,2,1);
idx = strcmp(cluster_infos.regionname, 'AM');
nsu = sum(strcmp('SU', cluster_infos{idx, 'clustype'}));
PlotTuning(tc, idx, xl, xt, yl, yt, alphacrit)
title(sprintf('%s (n=%d, %d SU)', 'AM', sum(idx), nsu));
ylabel('^{FR}/_{FRmax}');
xlabel('stim rank');

subplot(1,2,2);
idx = strcmp(cluster_infos.regionname, 'HC') | ...
    strcmp(cluster_infos.regionname, 'EC') | ...
    strcmp(cluster_infos.regionname, 'PHC');
nsu = sum(strcmp('SU', cluster_infos{idx, 'clustype'}));
PlotTuning(tc, idx, xl, xt, yl, yt, alphacrit)
title(sprintf('%s (n=%d, %d SU)', 'HC, EC, PHC', sum(idx), nsu));
ylabel('');
xlabel('');

print(figh, ...
    sprintf('Figure2BTuningCurvesPerMin%dResps.png', ...
    min_n_re_stim),'-dpng', '-r600')
