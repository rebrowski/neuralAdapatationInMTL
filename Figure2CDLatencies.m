clearvars -except sessions eeg datadir stimdir% do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

load('priming_latencies.mat');

regions(1).name = 'AM';
regions(2).name = 'EC';
regions(3).name = 'HC';
regions(4).name = 'PHC';
nregions = numel(regions);

%% plot latencies
for fi = 1:2
    
    clear idx
    
    if fi == 1
        measure = aggregatelatencies;
        keepidx = 1:size(aggregatelatencies,1);
        figprefix = 'Figure2D_Latency';
        ylabeltext = 'Latency (ms)';
    else
        measure = aggregatedurations;
        keepidx = [measure(:,1)>1 & measure(:,2) > 1]';
        figprefix = 'Figure2C_BurstDuration';
        ylabeltext = 'Burst-Duration (ms)';
    end

    figh = figure('color', 'w', 'visible', 'on');
    figh.PaperUnits = 'inches';
    figh.PaperPosition = [0 0 7.4 5.2];
    % display it somewhat similar to what will be plotted
    figh.Position = [200 200  figh.PaperPosition(3)*150 figh.PaperPosition(4)*150];
    fontSize = 8;
    fontSizeSmall = 6;


    nrows = 2;
    ncols = 3;

    for r=1:nregions + 1

        if r <= nregions 
            idx = find(strcmp(regions(r).name, latlookup.region) & keepidx);
            regname = regions(r).name;
        else
            idx = find(~strcmp('other', latlookup.region) & keepidx);
            regname = 'overall';
        end
        
        % plot diffs
        subplot(nrows, ncols, r)

        for ix = 1:numel(idx)
            hold on;
            plot(measure(idx(ix),:), '-k')
        end 
        xlim([0.5, 2.5]);
        ylim([0, 1000]);
        
        set(gca, 'XTick', [1 2]);
        set(gca, 'XTickLabel', {'Primed', 'Control'})
        set(gca, 'XTickLabelRotation', 45)
        
        
        ylabel(ylabeltext)
        title(sprintf('%s (n=%d)', regname, numel(idx)));
        
        % print stats
        
        [h pt ci tstats] = ttest(measure(idx,1), ...
                                 measure(idx,2));

        [psr statssr] = signrank(measure(idx,1), ...
                                 measure(idx,2));
        
        statstext1 = sprintf(['tt p vs. c M(SD) = ', ...
                            '%.2f(%2.f) vs. %.2f(%.2f); t= %.3g, p = %.3g'], ...
                             nanmean(measure(idx,1)), ...
                             nanstd(measure(idx,1)), ...
                             nanmean(measure(idx,2)), ...
                             nanstd(measure(idx,2)), ...
                             tstats.tstat, ...
                             pt);
        text(0.6, 900, textwrap({statstext1}, 40), 'FontSize', fontSizeSmall);

        
        statstext2 = sprintf(['sr p vs.c Md(IQR) = ', ...
                            '%.2f(%2.f) vs. %.2f(%.2f); p = %.3g'], ...
                             nanmedian(measure(idx,1)), ...
                             iqr(measure(idx,1)), ...
                             nanmedian(measure(idx,2)), ...
                             iqr(measure(idx,2)), ...
                             psr);
        text(0.6, 750, textwrap({statstext2}, 40), 'FontSize', fontSizeSmall);

        
    end

    % do a boxplot
    subplot(nrows, ncols, nrows*ncols)

    %rearrange data for iosr
    pidx = find(~strcmp('other', latlookup.region) & keepidx);

    lat = [measure(pidx,1); measure(pidx,2)];
    reg = [latlookup.region(pidx), latlookup.region(pidx)]';
    cond = [repmat(1, numel(measure(pidx,1)), 1); ...
            repmat(2, numel(measure(pidx,1)), 1)];


    [y,x,g] = iosr.statistics.tab2box(reg, lat,  cond);


    h = iosr.statistics.boxPlot(x,y,...
                                'scalewidth',false,'xseparator',false,...
                                'groupLabels',g,'showLegend',true, ...
                                'theme', 'colorlines', ... %'colorall' 'default'
                                'sampleSize', false, ...
                                'showLegend', false, ...
                                'notch', false, ...
                                'outlierSize', 4, ...
                                'symbolMarker', 'x',...
                                'showScatter', true,...
                                'scatterMarker', 'o', ...
                                'scatterSize', 4, ...
                                'showViolin', false, ...
                                'lineWidth', 0.3 ...
                                );

 
    ylabel(ylabeltext);
    if fi == 1
        ylim([0 800]);
    else
        ylim([0 500]);
    end
    print(figh, sprintf('%s%s%s_primed_vs_control.png', ...
                        'plots', filesep, figprefix), '-dpng', '-r600');

end
