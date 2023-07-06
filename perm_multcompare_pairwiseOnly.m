function [handles results colorbarhandles] = perm_multcompare(condA, condB,nperms, ...
                          clusteralpha, alpha, xaxis, plotcolors, ...
                                              conditionnames, ax1, ...
                                              markSize, LfontSiz)


%% function to run pairwise comparision and one-sample t-test
% against zero and plot all results into the same plot
%function [handles] = perm_multcompare(condA, condB, nperms, ...
%                          clusteralpha, ...
%                          alpha, xaxis, plotcolors,
%                          conditionnames, ax1, markSize, LfontSiz)


if ~exist('nperms', 'var')
    nperms = 1000;
end

if ~exist('clusteralpha', 'var')
    clusteralpha = 0.005;
end

if ~exist('alpha', 'var')
    alpha = 0.05;
end

if ~exist('xaxis', 'var')
    xaxis = 1:size(condA, 2);
end

if ~exist('plotcolors', 'var')
    plotcolors = [0.1059 0.6196 0.4667;
                   0.8510 0.3725 0.0078];
end

if ~exist('conditionnames', 'var')
    conditionnames = {'A', 'B'};
end

if ~exist('markSize', 'var')
    markSize = 10;
end

if ~exist('LfontSiz', 'var')
    LfontSiz = markSize;
end


results = '';
pmax = 0.2; % for p-value colorbar

%% run the paired t-test
doplot = false;
[pos_true neg_true pos_dist neg_dist mask hdls pvals_true] = ...
    perm_ttest(condA, condB, nperms, [],clusteralpha, alpha, doplot);

%% plot curves
if ~exist('ax1', 'var') | isempty(ax1)
    ax1 = gca;
end

lineh(1) = plot_signals(condA, xaxis, plotcolors(1,:));
lineh(2) = plot_signals(condB, xaxis, plotcolors(2,:));
xlim([xaxis(1) xaxis(end)]);



grid on
% make some space for things at the bottom
margin = 0.2;
pos = get(ax1, 'Position');

newheight = pos(4) * 0.7;
newspace = pos(4) * 0.3;
pos(2) = pos(2) + pos(4) - newheight; 
pos(4) = newheight;
set(ax1, 'Position', pos);



[lh objh] = legend(lineh, conditionnames, 'Location', ...
                   'North', 'Orientation', 'horizontal', 'FontSize', ...
                   LfontSiz); 
legend('boxoff');
set(objh,'linewidth',2);
% make sure the legend is not written on top of data
yl = ylim;
ylim([yl(1) yl(2)+diff(yl)*0.1]);

%% plot p-values of paired test
set(ax1, 'Color', 'none');
cbpos = [pos(1) pos(2)-0.86*newspace pos(3) pos(4)*0.1];
ax2 = axes('Position', cbpos);

try 
    colormap(ax2, flipud(parula));
catch err
    
    load('parula_colormap'); 
    prla = parula_colormap;
    colormap(ax2, flipud(prl));
end

psh = imagesc(xaxis, 1, pvals_true); caxis([0 pmax]); 
axis off; axis  xy;
top_ = cbpos(2) + cbpos(4);

%% annotate the test
legendpos1 = axes('Position', [cbpos(1)- 0.1,
                  cbpos(2),
                  0.1,
                  cbpos(4)]);
hold on

 plot(0.2,1, 'o', 'color', plotcolors(1,:), 'MarkerSize', markSize, ...
     'MarkerFaceColor', plotcolors(1,:));

 plot(0.8,1, 'o', 'color', plotcolors(2,:), 'MarkerSize', markSize, ...
     'MarkerFaceColor', plotcolors(2,:));
text(0.5, 1, 'vs.', 'HorizontalAlignment', 'center')
xlim([0 1]);
ylim([0.5 1.5]);
axis off;
%% annotate sign clusters paired test
anotpos = [cbpos(1), ...
           cbpos(2)+cbpos(4), ...
           cbpos(3), ...
           cbpos(4)*0.9];
ax3 = axes('Position', anotpos);
set(ax3, 'color', 'none')
xlim([xaxis(1) xaxis(end)]);
ylim([0 1])

axis off
for c = 1:length(pos_true.froms)
    k = patch(xaxis([pos_true.froms(c) pos_true.tos(c)...
                     pos_true.tos(c) pos_true.froms(c)]), ...
              [0 0 1 1], [0.7 0.7 0.7]);
    set(k, 'EdgeColor', 'none');
    x_ = xaxis(round(mean([pos_true.tos(c) pos_true.froms(c)])));
    text(x_,0.5, '*', 'HorizontalAlignment', 'center');
    if c == 1
        results = sprintf('%s > %s ', conditionnames{1}, conditionnames{2});
    end
    results = [results sprintf('from %d ms to %d ms, clustersize = %d, sum of t-values = %0.5g, p = %0.5g; ', ...
                               round(xaxis(pos_true.froms(c))), round(xaxis(pos_true.tos(c))), ...
                               pos_true.sizes(c), pos_true.sumts(c), pos_true.cluster_p(c))];
    
end

for c = 1:length(neg_true.froms)
    k = patch(xaxis([neg_true.froms(c) neg_true.tos(c)...
                     neg_true.tos(c) neg_true.froms(c)]), ...
              [0 0 1 1], [0.7 0.7 0.7]);
    set(k, 'EdgeColor', 'none');
    x_ = xaxis(round(mean([neg_true.tos(c) neg_true.froms(c)])));

    text(x_, 0.5, '*', 'HorizontalAlignment', 'center'); 
    if c == 1
        results = [results sprintf('%s < %s ', conditionnames{1}, conditionnames{2})];
    end
    results = [results sprintf('from %d ms to %d ms, clustersize = %d, sum of t-values = %0.5g, p = %0.5g; ', ...
                               round(xaxis(neg_true.froms(c))), round(xaxis(neg_true.tos(c))), ...
                               neg_true.sizes(c), neg_true.sumts(c), neg_true.cluster_p(c))];

end

if length(pos_true.froms) == 0 & length(neg_true.froms) == 0
    results = [results sprintf('%s vs. %s \\t n.s.\\n', conditionnames{1}, ...
               conditionnames{2})];
end


handles = [ax1 ax2 ax3 legendpos1 lh];

colorbarhandles = ax2;