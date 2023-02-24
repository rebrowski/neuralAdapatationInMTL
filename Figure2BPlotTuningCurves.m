clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

min_n_re_stim = 4; % inclusion criterion for units responding to at least min_n_re_stimuli
load(sprintf('tuningCurvesMin%dResponsesPerUnit.mat', min_n_re_stim));
disp(sprintf('Units with at least %d significant responses', min_n_re_stim));
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
disp('AM')
PlotTuning(tc, idx, xl, xt, yl, yt, alphacrit)

% save things into source data file for ncomms
sourceDataFileName = ['source_data_files_ncomms', filesep, 'SourceDataFigure2B.xlsx'];
t1.data = [tc(idx,1:xl(2), 1);tc(idx,1:xl(2),2)];
t1.condlabels = repmat({'primed'}, sum(idx), 1);
t1.condlabels = [t1.condlabels;repmat({'control'}, sum(idx), 1)];
t1.cellId = [1:sum(idx), 1:sum(idx)]';
writetable(struct2table(t1), sourceDataFileName, 'Sheet', 'AM');

% annotate the plot
title(sprintf('%s (n=%d, %d SU)', 'AM', sum(idx), nsu));
ylabel('^{FR}/_{FRmax}');
xlabel('stim rank');

subplot(1,2,2);
idx = strcmp(cluster_infos.regionname, 'HC') | ...
    strcmp(cluster_infos.regionname, 'EC') | ...
    strcmp(cluster_infos.regionname, 'PHC');
nsu = sum(strcmp('SU', cluster_infos{idx, 'clustype'}));
disp('HC, EC, PHC')
PlotTuning(tc, idx, xl, xt, yl, yt, alphacrit)

t2.data = [tc(idx,1:xl(2), 1);tc(idx,1:xl(2),2)];
t2.condlabels = repmat({'primed'}, sum(idx), 1);
t2.condlabels = [t2.condlabels;repmat({'control'}, sum(idx), 1)];
t2.cellId = [1:sum(idx), 1:sum(idx)]';
writetable(struct2table(t2), sourceDataFileName, 'Sheet', 'HC-EC-PHC');


title(sprintf('%s (n=%d, %d SU)', 'HC, EC, PHC', sum(idx), nsu));
ylabel('');
xlabel('');

print(figh, ...
    sprintf('%s%sFigure2BTuningCurvesPerMin%dResps.png', ...
    'plots', filesep, min_n_re_stim),'-dpng', '-r600')
