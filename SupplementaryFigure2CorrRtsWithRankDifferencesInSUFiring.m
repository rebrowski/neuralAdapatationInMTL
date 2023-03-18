

%% this script correlates differences in firing to second- and thirdmost response-eliciting stimuli to ...
%% reaction time differences of primed vs. control condition
clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

load('tuningCurvesMin1ResponsesPerUnitSfig2.mat');
load('reactiontimes_primed_control_category.mat');


% get the average diff of primed - control of response magnitude (normalized to most-response-elicitng stim in the control condition)...
% at rank 2 for each subject

usubs = unique(subid_wide);
ususb2 = unique(cluster_infos.subjid);

regions = {'AM'; 'otherMTL'};

h = figure('Color', 'w', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 8]);

for r = 1:numel(regions)

    if strcmp(regions{r}, 'AM')
        regIdx = strcmp(cluster_infos.regionname, 'AM');
    elseif strcmp(regions{r}, 'otherMTL')
        regIdx = strcmp(cluster_infos.regionname, 'EC') | ...
            strcmp(cluster_infos.regionname, 'PHC') | ...
            strcmp(cluster_infos.regionname, 'HC');
    end

    for s = 1:numel(usubs)

        % get average RT diff across sessions for a subject
        idx1 = (subid_wide == usubs(s))';
        rtprimed(s) = nanmean(nanmean(rtcat_wide_sr(idx1,:,1)));
        rtcontrol(s) = nanmean(nanmean(rtcat_wide_sr(idx1,:,2)));
        rtdiff(s) = rtprimed(s) - rtcontrol(s);

        % get average diff of unit activity across units for a subject
        idx2 = cluster_infos.subjid == usubs(s) & regIdx;

        rank1diff(s) = mean(tc(idx2,1,1)) - mean(tc(idx2,1,2));
        rank2diff(s) = mean(tc(idx2,2,1)) - mean(tc(idx2,2,2));
        rank3diff(s) = mean(tc(idx2,3,1)) - mean(tc(idx2,3,2));

    end

    % for one subject, there are no units with at least four responses in that region.
    idBoth = ~isnan(rank1diff);

    subplot(2,2,(r-1)*2 + 1);
    scatter(rtdiff(idBoth), rank2diff(idBoth), 'filled');
    xlabel('RT(primed-control, seconds)')
    ylabel({'Normalized Firing Diff', 'Rank 2 (primed - control)'});
    [rho p] = corr(rtdiff(idBoth)', rank2diff(idBoth)', 'Type', 'Spearman');
    title(sprintf('%s, Spearman R = %.2f, p = %.2g', regions{r}, rho, p ),'FontWeight', 'normal');

    
    if r == 1
        t.rtdiffAmy = rtdiff(idBoth)';
        t.rank2DiffAmy = rank2diff(idBoth)';
    elseif r == 2
        t.rtdiffOther = rtdiff(idBoth)';
        t.rank2DiffOther = rank2diff(idBoth)';
    end

    subplot(2,2,(r-1)*2 +2);
    scatter(rtdiff(idBoth), rank3diff(idBoth), 'filled');
    xlabel('RT(primed-control, seconds)')
    ylabel({'Normalized Firing Diff', 'Rank 3 (primed - control)'});
    [rho p] = corr(rtdiff(idBoth)', rank3diff(idBoth)', 'Type', 'Spearman');
    title(sprintf('%s, Spearman R = %.2f, p = %.2g', regions{r}, rho, p ), 'FontWeight', 'normal');


    if r == 1
        t.rank3DiffAmy = rank2diff(idBoth)';
    elseif r == 2
        t.rank3DiffOther = rank2diff(idBoth)';
    end

end

save('rankdiffsPerSubject.mat', 'rank1diff', 'rank2diff', 'rank3diff'); % to be used in ospr_repetitions_reviewer1_2.m that also looks at correlations of SU Firing with iEEG activity
sourceDataFileName = ['source_data_files_ncomms', filesep ,'SourceDataSupplementaryFigure2.xlsx'];

t.rtdiffAmy(end+1) = NaN;
t.rank2DiffAmy(end+1) = NaN;
t.rank3DiffAmy(end+1) = NaN;
writetable(struct2table(t), sourceDataFileName);
print(h, ['plots',filesep, 'SupplementaryFigure2.png'], '-dpng');

