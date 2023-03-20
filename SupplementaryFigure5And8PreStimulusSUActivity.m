% assess evidence for spreading activation, first calculate zscores also in
% the baseline separate per condition, do a t-test in each region whether
% activity is higher for primed vs. control stimuli


clearvars -except sessions

z = load('zvals_condition.mat'); % computed by ospr_calculate_zscores_condition.m
cr = load('category_responses.mat'); % get responsive units (priming & control collapsed), user consider_rs (ranksum) to be consistent with 
min_n_re_stim = 1; % take all units with at least one response
uidxcombined = sum(cr.consider_rs,2) >= min_n_re_stim; % p/c combined


regions(1).name='AM';
regions(2).name='otherMTL';

% initialize a struct with aggregated data

t.region = {};
t.condition = {};
t.zscorePreStim = [];
t.zscorePostStim = [];
t.frPreStim = [];
t.stimName = {};

cc = 0

for r = 1:numel(regions)
    % find units in region
    if strcmp(regions(r).name, 'AM')
        uidxInRegion = uidxcombined & strcmp(cr.cluster_lookup.regionname, 'AM');

    elseif strcmp(regions(r).name, 'otherMTL')
        uidxInRegion = uidxcombined & (strcmp(z.cluster_lookup.regionname, 'HC') | ...
            strcmp(z.cluster_lookup.regionname, 'EC') | ...
            strcmp(z.cluster_lookup.regionname, 'PHC'));
    end
    
    fuidx = find(uidxInRegion);
    
    for u = 1:numel(fuidx)
        % get most response-eliciting stim
        [m_ stimId] = min(cr.pvals_rs(fuidx(u),:));
        stimName = cr.stim_lookup(stimId);

        for c = 1:2
            cc = cc+1;
            t.region = [t.region; regions(r).name];
            t.stimName = [t.stimName; stimName];
            if c == 1
                t.condition = [t.condition; 'primed'];
            elseif c == 2
                t.condition = [t.condition; 'control'];
            end
            t.zscorePostStim = [t.zscorePostStim; z.zvals(fuidx(u), stimId, c)];
            t.zscorePreStim = [t.zscorePreStim; z.zvalsPreStim(fuidx(u), stimId, c)];
            t.frPreStim = [t.frPreStim; z.frPreStim(fuidx(u), stimId, c)];
        end
    end
end

%% do a figure with raw firing rates


h1 = figure('Color', 'w', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 5]);
subplot(1,2,1);
r1idx = strcmp(t.region, 'AM');
boxplot(t.frPreStim(r1idx), t.condition(r1idx));
box off
ylabel('Pre-Stimulus FR (Hz)')
title(sprintf('AM n=%d', sum(r1idx)/2)); % divide indices by two because we have two entries for each neuron

primedIdx = strcmp(t.condition, 'primed');
controlIdx = strcmp(t.condition, 'control');

[h p ci stats] = ttest(t.frPreStim(r1idx & primedIdx), t.frPreStim(r1idx & controlIdx));
text(1,36, sprintf('t_{%d}=%.3g, p=%.3g', stats.df, stats.tstat, p));
disp('AM')
disp('========')
disp(sprintf('M(SD)_{primed} = %.2g(%.2g), M(SD)_{control} = %.2g(%.2g)', ...
    mean(t.frPreStim(r1idx & primedIdx)), std(t.frPreStim(r1idx & primedIdx)), ...
    mean(t.frPreStim(r1idx & controlIdx)), std(t.frPreStim(r1idx & controlIdx))));
disp(sprintf('signrank p = %.3g', signrank(t.frPreStim(r1idx & primedIdx), t.frPreStim(r1idx & controlIdx))));

subplot(1,2,2);
r2idx = strcmp(t.region, 'otherMTL');
boxplot(t.frPreStim(r2idx), t.condition(r2idx));
box off
ylabel('Pre-Stimulus FR (Hz)')
title(sprintf('otherMTL n=%d', sum(r2idx)/2));
[h p ci stats] = ttest(t.frPreStim(r2idx & primedIdx), t.frPreStim(r2idx & controlIdx));

text(1,36, sprintf('t_{%d}=%.3g, p=%.3g', stats.df, stats.tstat, p));
disp('otherMTL')
disp('========')
disp(sprintf('M(SD)_{primed} = %.2g(%.2g), M(SD)_{control} = %.2g(%.2g)', ...
    mean(t.frPreStim(r2idx & primedIdx)), std(t.frPreStim(r2idx & primedIdx)), ...
    mean(t.frPreStim(r2idx & controlIdx)), std(t.frPreStim(r2idx & controlIdx))));
disp(sprintf('signrank p = %.3g', signrank(t.frPreStim(r2idx & primedIdx), t.frPreStim(r2idx & primedIdx))));
print(h1, ['plots' filesep  'SupplementaryFigure8Hz.png'], '-dpng');


%% do a figure for Z-scores 

h2 = figure('Color', 'w', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 5]);
subplot(1,2,1);
r1idx = strcmp(t.region, 'AM');
boxplot(t.zscorePreStim(r1idx), t.condition(r1idx));
box off
ylabel('Z(FR)')
title(sprintf('AM n=%d', sum(r1idx)/2)); % divide indices by two because we have two entries for each neuron

primedIdx = strcmp(t.condition, 'primed');
controlIdx = strcmp(t.condition, 'control');

[h p ci stats] = ttest(t.zscorePreStim(r1idx & primedIdx), t.zscorePreStim(r1idx & controlIdx));
text(1,2, sprintf('t_{%d}=%.3g, p=%.3g', stats.df, stats.tstat, p));
disp('AM')
disp('========')
disp(sprintf('M(SD)_{primed} = %.2g(%.2g), M(SD)_{control} = %.2g(%.2g)', ...
    mean(t.zscorePreStim(r1idx & primedIdx)), std(t.zscorePreStim(r1idx & primedIdx)), ...
    mean(t.zscorePreStim(r1idx & controlIdx)), std(t.zscorePreStim(r1idx & controlIdx))));



subplot(1,2,2);
r2idx = strcmp(t.region, 'otherMTL');
boxplot(t.zscorePreStim(r2idx), t.condition(r2idx));
box off
ylabel('Z(FR)')
title(sprintf('otherMTL n=%d', sum(r2idx)/2));
[h p ci stats] = ttest(t.zscorePreStim(r2idx & primedIdx), t.zscorePreStim(r2idx & controlIdx));

text(1,2, sprintf('t_{%d}=%.3g, p=%.3g', stats.df, stats.tstat, p));
disp('otherMTL')
disp('========')
disp(sprintf('M(SD)_{primed} = %.2g(%.2g), M(SD)_{control} = %.2g(%.2g)', ...
    mean(t.zscorePreStim(r2idx & primedIdx)), std(t.zscorePreStim(r2idx & primedIdx)), ...
    mean(t.zscorePreStim(r2idx & controlIdx)), std(t.zscorePreStim(r2idx & controlIdx))));

print(h2, ['plots' filesep 'SupplementaryFigure5Zscores.png'], '-dpng');

writetable(struct2table(t), ['source_data_files_ncomms' filesep 'SourceDataSupplementaryFigure5and8.xlsx']);

