%% aggregate data to calculate an ANOVA of Region X Rank (1 vs. 2) X Condition (primed vs. control)

clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end
load('tuningCurvesMin4ResponsesPerUnit.mat');
regions = {'AM'; 'otherMTL'};

% prepare output

t.fr = [];
t.zfr = [];
t.rank = [];
t.condition = {};
t.region = {};
t.subject = [];
t.sessid = [];
t.unitId = [];

oc = 1;

for r = 1:numel(regions)
    
    if strcmp(regions{r}, 'AM')
        regIdx = strcmp(cluster_infos.regionname, 'AM');
        regname = 'AM';
        
    elseif strcmp(regions{r}, 'otherMTL')
        regIdx = strcmp(cluster_infos.regionname, 'EC') | ...
            strcmp(cluster_infos.regionname, 'PHC') | ...
            strcmp(cluster_infos.regionname, 'HC');
        regname = 'otherMTL';
    end

    UnitsinRegion = find(regIdx);

    for ui = 1:numel(UnitsinRegion)
        for rank = 1:4
            for condition = 1:2
                cond_ = 'primed'; 
                if condition == 2; cond_ = 'control'; end;
                t.fr = [t.fr; tcfr(UnitsinRegion(ui),rank,condition)];
                t.zfr = [t.zfr; tc(UnitsinRegion(ui),rank,condition)];
                t.rank = [t.rank; rank];
                t.condition = [t.condition; cond_];
                t.region = [t.region; regname];
                t.subject = [t.subject; cluster_infos.subjid(UnitsinRegion(ui))];
                t.sessid = [t.sessid; cluster_infos.sessid(UnitsinRegion(ui))];
                t.unitId = [t.unitId; UnitsinRegion(ui)];
            end
        end
    end
end

tb = struct2table(t);
writetable(tb, 'FiringRatesByRankCondition.csv');


