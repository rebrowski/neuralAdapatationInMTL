function PlotTuning(tc, idx, xl, xt, yl, yt, alphacrit)
%% helper function for ospr_plot_tuning_curves

mrksz = 15;% 'MarkerSize'


plot_signals(tc(idx,1:xl(2), 1), 1:xl(2), 'b');
plot(1:xl(2), mean(tc(idx,1:xl(2), 1)), '.b', 'MarkerSize', mrksz);

plot_signals(tc(idx,1:xl(2), 2), 1:xl(2), 'r');
plot(1:xl(2), mean(tc(idx,1:xl(2), 2)), '.r', 'MarkerSize', mrksz);

%% anonate stats
prs = NaN(1, xl(2));
for x = 1:xl(2)
    
    [k_ prs(x) CI stats] = ttest(tc(idx,x,1), tc(idx,x,2));
    
    disp(sprintf('rank %d primed vs. control: t(%d)=%.2g, CI=[%.2g %.2g], p = %.2g', x, stats.df, stats.tstat, ...
        CI(1), CI(2), prs(x)));

    fw = 'normal';
    if prs(x) < alphacrit
        issigstr = '*';
        fw = 'bold';
    end
    
    text(x, yl(1), sprintf('p=%.2g%s', prs(x)), ...
        'Rotation', 90, 'VerticalAlignment', 'middle', ...
        'HorizontalAlignment', 'left', 'FontWeight', fw);
end

%% calculate t-tests b
diffs = tc(idx,1:xl(2), 1) - tc(idx, 1:xl(2), 2);
for id = 1:xl(2)-1
    diffdiffs(id, :) = diffs(:,id+1) - diffs(:,id);
    
    [h prs(id) CI stats] = ttest(diffs(:,id), diffs(:,id+1));
    
    fw = 'normal';
    if prs(id) < alphacrit
        issigstr = '*';
        fw = 'bold';
    end
    
    %disp(sprintf('rank %d vs. %d: t(%d)=%.2g, CI=[%.2g %.2g], p = %.2g', id, id+1, stats.df, stats.tstat, ...
       % CI(1), CI(2), prs(id)));
end

grid on
set(gca, 'xlim', xl);
set(gca, 'ylim', yl);
set(gca, 'XTick', xt);
set(gca, 'YTick', yt);



