%% this script correlations reaction time differences of primed vs. control condition with ERP differences at different time points
clearvars -except sessions eeg datadir stimdir
if ~exist('datadir', 'var')
    startup % the one in this repo
end

% params: 
tw1from = 200; tw1to = 350;
tw2from = 350; tw2to = 600;

% load behavior data
rts = load('reactiontimes_primed_control_category.mat');

% load erps
if ~exist('eeg', 'var')
    load([datadir filesep 'ieeg_linked_mastoids_256Hz.mat']);
end

%% go through ERPS an calculate mean amplitudes per subject in two
% time-windows
nsessions = numel(eeg);

% load regions info an add some info for ieeg
regions(1).ieegsites = {'AL', 'AR'};
regions(1).sites = {'LA', 'RA'};
regions(1).name = 'AM';

regions(2).sites = {'LAH', 'RAH', 'LMH', 'RMH', 'LPH', 'RPH', 'LEC', 'REC', 'LPHC', 'RPHC'};
regions(2).ieegsites = {'AHL', 'AHR', 'MHL', 'MHR', 'PHL', 'PHR', 'ECL', 'ECR', 'PHCL', 'PHCR'};
regions(2).name = 'otherMTL';

% aggregate data  

for sessi = 1:nsessions
    
    % get indices to samples for timewindows
    tw1idx = eeg(sessi).stime > tw1from & eeg(sessi).stime <= tw1to;
    tw2idx = eeg(sessi).stime > tw2from & eeg(sessi).stime <= tw2to;

    for regi = 1:numel(regions)

        % 1. get channels matching region
        nsites_in_region = numel(regions(regi).ieegsites);
        siteidx = false(1,eeg(sessi).nsites);
        for sitei = 1:nsites_in_region
            i_ = find(strcmp(regions(regi).sites(sitei), ...
                eeg(sessi).sites));
            if ~isempty(i_)
                siteidx(i_) = true;
            end

        end

        if sum(siteidx) > 0
            siteidxf = find(siteidx);
            clear dat

            for chani = 1:numel(siteidxf)

                for cond = 1:2

                    tidx = eeg(sessi).condition.condition == cond & ...
                        ~eeg(sessi).isartefact(siteidxf(chani),:);

                    % 1. average over trials per channel
                    dat(chani, cond, :) ...
                        = squeeze(mean(eeg(sessi).sdata(siteidxf(chani),tidx,:),2));

                end % conditions
            end % channels in site

            % 2. average over channels, and timewindows defined above

            for cond = 1:2
                % erp
                erp(regi).dat(sessi,cond, :) = ...
                    squeeze(mean(dat(:,cond,:),1));

                % timewindows                
                erp(regi).tw1(sessi, cond) = mean(erp(regi).dat(sessi, cond,tw1idx));
                erp(regi).tw2(sessi, cond) = mean(erp(regi).dat(sessi, cond,tw2idx));

            end

            erp(regi).subject(sessi) = eeg(sessi).condition.subject;
        end % if


    end % regions
    erp(regi).regionname = regions(regi).name;

    
end % sessions

% plot ERPs for each region (sanity check)

h1 = figure('Color', 'w', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 8]);

for regi = 1:2

    subplot(2,1,regi)
    plot(eeg(1).stime, squeeze(mean(erp(regi).dat(:, 1,:), 1))); hold on;
    plot(eeg(1).stime, squeeze(mean(erp(regi).dat(:, 2,:), 1))); hold on;
    
end

% average over subjects, compute difference primed - control, and do a
% figure
assert(sum(erp(regi).subject - rts.subid_wide) == 0); % make sure, RT & ERP data corresponds
usubs = unique(rts.subid_wide);


for r = 1:numel(regions)

    for subji = 1:numel(usubs)

        % RT primed - control
        idx1 = (rts.subid_wide == usubs(subji))';
        rtprimed(subji) = nanmean(nanmean(rts.rtcat_wide_sr(idx1,:,1)));
        rtcontrol(subji) = nanmean(nanmean(rts.rtcat_wide_sr(idx1,:,2)));
        rtdiff(subji) = rtprimed(subji) - rtcontrol(subji);

        subidx = erp(r).subject == usubs(subji);

        % ERP primed - control
        tw1(r).avErp(subji) = mean(erp(r).tw1(subidx, 1))-mean(erp(r).tw1(subidx, 2));
        tw2(r).avErp(subji) = mean(erp(r).tw2(subidx, 1))-mean(erp(r).tw2(subidx, 2));

    end

end

%% plot correlations of behavior with iEEG

h = figure('Color', 'w', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 8]);
for r = 1:2

    % plot correlations with first and second timewindow
    subplot(2,2,(r-1)*2 + 1);
    scatter(rtdiff, tw1(r).avErp, 'filled');
    xlabel('RT(primed-control, seconds)')
    ylabel({'iEEG ERP(primed-control), mV', sprintf('average over %d to %d ms', tw1from, tw1to)});
    [rho p] = corr(rtdiff', tw1(r).avErp', 'Type', 'Spearman');
    title(sprintf('%s, Spearman R = %.2f, p = %.2g', regions(r).name, rho, p ),'FontWeight', 'normal');


    subplot(2,2,(r-1)*2 +2);
    scatter(rtdiff, tw2(r).avErp, 'filled');
    xlabel('RT(primed-control, seconds)')
    ylabel({'iEEG ERP(primed-control), mV', sprintf('average over %d to %d ms', tw2from, tw2to)});
    [rho p] = corr(rtdiff', tw2(r).avErp', 'Type', 'Spearman');
    title(sprintf('%s, Spearman R = %.2f, p = %.2g', regions(r).name, rho, p ),'FontWeight', 'normal');

    
    if r == 1
        sourcedatatable.rtdiffAmy = rtdiff';
        sourcedatatable.tw1AmyiEEGDiff = tw1(r).avErp';
        sourcedatatable.tw2AmyiEEGDiff = tw2(r).avErp';
    else
        sourcedatatable.rtdiffOther = rtdiff';
        sourcedatatable.tw1OtheriEEGDiff = tw1(r).avErp';
        sourcedatatable.tw2OtheriEEGDiff = tw2(r).avErp';
    end

end


print(h, ['plots',filesep, 'SupplementaryFigure3CorrRTsiEEG.png'], '-dpng');

writetable(struct2table(sourcedatatable), ['source_data_files_ncomms' filesep 'SourceDataSupplementaryFigure3.xlsx']);


% last point of Reviewer 1: correlate unit firing with iEEG, so lets load
% what we did in response to the first point of reviewer 1, i.e.,
% SupplementaryFigure2CorrRtsWithRankDiffsinSUFiring.m
su = load('rankdiffsPerSubject.mat');

subidx = ~isnan(su.rank1diff); % there is one subject having NaNs

h2 = figure('Color', 'w', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 8]);
for r = 1:2

    if r == 1 % amygdala: Sharpening rank 2 diffs
        dv = su.rank2diff(subidx);
        label = {'SU Activity difference', 'primed-control for rank 2 stimuli'};
    elseif r == 2 % other regions: fatiguing rank 1 diffs
        dv = su.rank1diff(subidx);
        label = {'SU Activity difference', 'primed-control for rank 1 stimuli'};
    end

    % plot correlations with first and second timewindow
    subplot(2,2,(r-1)*2 + 1);
    scatter(dv, tw1(r).avErp(subidx), 'filled');
    xlabel(label);
    ylabel({'iEEG ERP(primed-control), mV', sprintf('average over %d to %d ms', tw1from, tw1to)});
    [rho p] = corr(dv', tw1(r).avErp(subidx)', 'Type', 'Spearman');
    title(sprintf('%s, Spearman R = %.2f, p = %.2g', regions(r).name, rho, p ),'FontWeight', 'normal');


    subplot(2,2,(r-1)*2 +2);
    scatter(dv, tw2(r).avErp(subidx), 'filled');
    xlabel(label);
    ylabel({'iEEG ERP(primed-control), mV', sprintf('average over %d to %d ms', tw2from, tw2to)});
    [rho p] = corr(dv', tw2(r).avErp(subidx)', 'Type', 'Spearman');
    title(sprintf('%s, Spearman R = %.2f, p = %.2g', regions(r).name, rho, p ),'FontWeight', 'normal');


    if r == 1
        sourcedatatable2.AmyRank2Diff=dv';
        sourcedatatable2.AmyTw1ERPs=tw1(r).avErp(subidx)';
        sourcedatatable2.AmyTw2ERPs=tw2(r).avErp(subidx)';
    else
        sourcedatatable2.OtherRank1Diff=dv';
        sourcedatatable2.OtherTw1ERPs=tw1(r).avErp(subidx)';
        sourcedatatable2.OtherTw2ERPs=tw2(r).avErp(subidx)';
    end
end

writetable(struct2table(sourcedatatable2), ['source_data_files_ncomms' filesep 'SourceDataSupplementaryFigure6.xlsx']);
print(h2, ['plots', filesep 'SupplementaryFigure6.png'], '-dpng');

