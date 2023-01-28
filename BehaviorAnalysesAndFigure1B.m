clearvars -except sessions

secondleveldir = '~/projects/ospr/secondlevel/';
% load behavior data, vars with _sr at the end are restricted to trials in
% which the prime in the control condtiion afforded the same response as
% the target (ie., was of a different category, but same meta-category)
load([secondleveldir, 'reactiontimes_primed_control_category.mat']);

%% 1. calcuate average over categories per session

t.rt_sr = []; % reaction times same response
t.rt_all = []; % reaction times overall
t.pc_sr = []; % percept correct same response
t.pc_all = []; % percent correct overall
t.cond = {}; % primed or control
t.subject = [];

for s = 1:numel(subid_wide)
    
    % primed condition
    t.cond = [t.cond; {'primed'}];
    t.subject = [t.subject; subid_wide(s)];    
    t.rt_all = [t.rt_all; median(rtcat_wide(s, :, 1), 'omitnan')*1000];
    t.rt_sr = [t.rt_sr; median(rtcat_wide_sr(s, :, 1), 'omitnan')*1000];
    t.pc_all = [t.pc_all; mean(pcorrcat_wide(s, :, 1), 'omitnan')*100];
    t.pc_sr = [t.pc_sr; mean(pcorrcat_wide_sr(s, :, 1), 'omitnan')*100];
    
    % control condition
    t.cond = [t.cond; {'control'}];    
    t.subject = [t.subject; subid_wide(s)];
    t.rt_all = [t.rt_all; median(rtcat_wide(s, :, 2),  'omitnan')*1000];
    t.rt_sr = [t.rt_sr; median(rtcat_wide_sr(s, :, 2),  'omitnan')*1000];
    t.pc_all = [t.pc_all; median(pcorrcat_wide(s, :, 2), 'omitnan')*100];
    t.pc_sr = [t.pc_sr; median(pcorrcat_wide_sr(s, :, 2), 'omitnan')*100];

end

%% print stats to console
primedIdx = strcmp(t.cond, 'primed');
controlIdx = strcmp(t.cond, 'control');

[P, H] = signrank(t.rt_all(primedIdx), t.rt_all(controlIdx)); 
fprintf('Med(IQR) RT primed overall:%.2f (%.2f), control %.2f (%.2f); signrank p = %.3g\n', ...
             median(t.rt_all(primedIdx)), iqr(t.rt_all(primedIdx)), median(t.rt_all(controlIdx)), iqr(t.rt_all(controlIdx)), P );
       
[P, H] = signrank(t.rt_sr(primedIdx), t.rt_sr(controlIdx)); 
fprintf('Med(IQR) RT primed same response:%.2f (%.2f), control %.2f (%.2f); signrank p = %.3g\n', ...
             median(t.rt_sr(primedIdx)), iqr(t.rt_sr(primedIdx)), median(t.rt_sr(controlIdx)), iqr(t.rt_sr(controlIdx)), P );

[P, H] = signrank(t.pc_all(primedIdx), t.pc_all(controlIdx)); 
fprintf('Med(IQR) percent correct primed overall:%.2f (%.2f), control %.2f (%.2f); signrank p = %.3g\n', ...
             median(t.pc_all(primedIdx)), iqr(t.pc_all(primedIdx)), median(t.pc_all(controlIdx)), iqr(t.pc_all(controlIdx)), P );
       
[P, H] = signrank(t.pc_sr(primedIdx), t.pc_sr(controlIdx)); 
fprintf('Med(IQR) percent correct primed same response:%.2f (%.2f), control %.2f (%.2f); signrank p = %.3g\n', ...
             median(t.pc_sr(primedIdx)), iqr(t.pc_sr(primedIdx)), median(t.pc_sr(controlIdx)), iqr(t.pc_sr(controlIdx)), P );



%% do a figure for RTs
h1 = figure('Color', 'w', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 4]);
ax1 = subplot(1,2,1);

boxplot(t.rt_all, t.cond);
title('all trials');
ylabel('RT(ms)');
box off

[P, H] = signrank(t.rt_all(primedIdx), t.rt_all(controlIdx)); 
%fprintf('Med(IQR) RT primed overall:%.2f (%.2f), control %.2f (%.2f); signrank p = %.3g\n', ...
%             median(t.rt_all(primedIdx)), iqr(t.rt_all(primedIdx)), median(t.rt_all(controlIdx)), iqr(t.rt_all(controlIdx)), P );
text(1.2,1400, sprintf('p=%.3g*', P));
hold on;
plot([1,2],[1300, 1300], '-k');

ax2 = subplot(1,2,2);
boxplot(t.rt_sr, t.cond);
title({'control targets primed', 'by same meta-category'});
ylabel('RT(ms)')
box off

[P, H] = signrank(t.rt_sr(primedIdx), t.rt_sr(controlIdx)); 
%fprintf('Med(IQR) RT primed same response:%.2f (%.2f), control %.2f (%.2f); signrank p = %.3g\n', ...
%             median(t.rt_sr(primedIdx)), iqr(t.rt_sr(primedIdx)), median(t.rt_sr(controlIdx)), iqr(t.rt_sr(controlIdx)), P );
text(1.2,1400, sprintf('p=%.3g*', P));
hold on;
plot([1,2],[1300, 1300], '-k');

linkaxes([ax1, ax2], 'y');

print(h1, 'Figure1B.png', '-dpng');