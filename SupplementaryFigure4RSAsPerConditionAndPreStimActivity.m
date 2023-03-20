% calculate neural similarity between two stimuli, look at tuning curves
% for more similar vs. dissimilar primes & targets
%% 

clearvars -except sessions ztrial
startup
% load stuff
if ~exist('sessions', 'var')
    load([datadir, 'sessions.mat']);
end

if ~exist('ztrial', 'var')
    ztrial = load([datadir filesep 'zvals_trials_and_pre_trial.mat']); % done in ospr_calculate_zscores_trials_and_prevtrials.m
end
z = load('zvals_100_1000.mat');
cr = load('category_responses.mat'); % get responsive units (priming & control collapsed), user consider_rs (ranksum) to be consistent with 

min_n_re_stim = 1; % take all units with at least one response
uidxcombined = sum(cr.consider_rs,2) >= min_n_re_stim; % p/c combined

regions(1).name='AM';
regions(2).name='otherMTL';

%% 1. calculate RDMs for all stimulus paris, save median of distance values
% prepare output for boxplot
distances.dist = [];
distances.region = [];
distances.condition = [];
distances.conditionInitial = [];

% get indices in stimulus RDM
% calculate median similarity for within vs. across category pairs
idxStimsInRdm = zeros(100,100);

for c1 = 1:100
    for c2 = c1:100
        if c1~= c2 && floor((c1-1)/10) == floor((c2-1)/10)
            idxStimsInRdm(c1,c2) = 2;
        elseif c1~= c2 && floor((c1-1)/10) ~= floor((c2-1)/10)
            idxStimsInRdm(c1,c2) = 1;
        end
    end
end

% look at similarities for primed, control condition overall in regions we
% are intereted in and make a figure for the reviewer
for r = 1:numel(regions)
    
    if strcmp(regions(r).name, 'otherMTL')
        idx = strcmp(z.cluster_lookup.regionname, 'HC') | ...
            strcmp(z.cluster_lookup.regionname, 'EC') | ...
            strcmp(z.cluster_lookup.regionname, 'PHC');

    elseif(strcmp(regions(r).name, 'AM'))
        idx = strcmp(z.cluster_lookup.regionname, 'AM');
    end

    S=pdist(z.zvals(idx, :)', 'correlation');    
    regions(r).RDM = squareform(S); % store RDM for further analyses below
    
    % store values for plotting
    distprimed = regions(r).RDM(idxStimsInRdm == 2);
    distprimedmedianstplit = distprimed > median(distprimed); % higher distance than medin == 1
    primedLowLabels = repmat({'low similarity primes'}, sum(distprimedmedianstplit == 1), 1);
    primedHighLabels = repmat({'high similarity primes'},sum(distprimedmedianstplit == 0),1);
    
    distcontrol = regions(r).RDM(idxStimsInRdm == 1);
    distcontrol = distcontrol(:);
    controllabels = repmat({'control'},numel(distcontrol),1);
    
    distances.dist = [distances.dist; distprimed(distprimedmedianstplit==0); distprimed(distprimedmedianstplit == 1); distcontrol];
    distances.condition = [distances.condition; primedHighLabels; primedLowLabels; controllabels];
    distances.conditionInitial = [distances.conditionInitial; repmat({'primed'}, numel(primedHighLabels)+numel(primedLowLabels),1); controllabels];
    distances.region = [distances.region; repmat({regions(r).name}, sum(sum(idxStimsInRdm > 0)),1)];

end


% 2. make plots for stimulus pair similarities by condition and region (see above
h1 = figure('Color', 'w', 'PaperUnits', 'inches', 'PaperPosition', [0 0 8 5]);
subplot(1,2,1);
r1idx = strcmp(distances.region, 'AM');
boxplot(distances.dist(r1idx), distances.conditionInitial(r1idx));
box off
ylabel('1-R')
title('AM')
ylim([0.5 1.1]);

subplot(1,2,2);
r2idx = strcmp(distances.region, 'otherMTL');
boxplot(distances.dist(r2idx), distances.conditionInitial(r2idx));
box off
ylabel('1-R')
title('otherMTL')
ylim([0.5 1.1]);

print(h1, ['plots' filesep 'SupplementaryFigure4StimSimilarityRSAPerCondition.png'], '-dpng');

writetable(struct2table(distances), ['source_data_files_ncomms' filesep 'SourceDataSupplementaryFigure4.xlsx']);

%% 2. look at repetition supression for prime-target pairs high in similarity and for the ones low in similarity
% the idea here is to go through each responsive unit, take the most
% response-eliciting stimulus and store z-scores during the response windwo
% for 1. control condition (5 trials), 2. high similarity primes, 3. low
% similarity primes

% define some vars for the output
t.unitNo = [];
t.zscores = [];
t.condition = {};
t.regionname = {};
t.prevStimName = {};
t.stimName = {};
t.similarityToPrevStim = []; %R-1 according to the above RDMs for the subregion we are in.
t.similarityCondition = {}; % group experimental (primed) and control (non-primed, diff category) in high and low similarity trials
%t.similarityRankWithinCondition = [];


% go through all units
nunits = sum(uidxcombined);

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

        % get trials idx for control, high and low similarity primed
        % condition
        tidx = strcmp(ztrial.stim_lookup, stimName);
        ctidx = tidx & ztrial.cond_lookup == 2;
        ptidx = tidx & ztrial.cond_lookup == 1; 

        % append data to output above
        t.unitNo = [t.unitNo; repmat(fuidx(u), sum(tidx), 1)]; % ten trials
        t.zscores = [t.zscores; ztrial.zvals(fuidx(u), ptidx)'; ztrial.zvals(fuidx(u), ctidx)'];
        t.condition = [t.condition; repmat({'primed'}, sum(ptidx), 1); repmat({'control'}, sum(ctidx), 1)];
        t.regionname = [t.regionname; repmat({regions(r).name}, sum(tidx), 1)];
        t.stimName = [t.stimName; repmat(stimName, sum(tidx), 1)];
        t.prevStimName = [t.prevStimName; ztrial.stimprevname(fuidx(u),ptidx)'; ztrial.stimprevname(fuidx(u), ctidx)'];
        
        % get R-1 for each of the ten trials above
        for tid = 1:10
            xi = find(strcmp(z.stim_lookup, stimName)); % col in RDM
            s_ = t.prevStimName{length(t.prevStimName)-10+tid};
            yi = find(strcmp(z.stim_lookup, s_(1:length(s_)-4))); % row in RDM
            t.similarityToPrevStim(length(t.prevStimName)-10+tid, 1) = regions(r).RDM(xi,yi);

        end
        
        % sort trials in exp and control and assign to high vs. low        
        % similarity prime target pairs

        % initialize entries
        t.similarityCondition = [t.similarityCondition; repmat({' '}, 10, 1)];

        % experimental condition
        [primeSimilarityP primeSimilarityIdxPrimed] = sort(t.similarityToPrevStim(end-9:end-5), 'ascend');
        % control condition
        [primeSimilarityC primeSimilarityIdxControl] = sort(t.similarityToPrevStim(end-4:end), 'ascend');
        for i_ = 1:5
            if i_ < 3                
                t.similarityCondition{end-10+primeSimilarityIdxPrimed(i_),1} = 'PrimedHighSimilarity';
                t.similarityCondition{end-5+primeSimilarityIdxControl(i_),1} = 'ControlHighSimilarity';
            elseif i_ > 3
                t.similarityCondition{end-10+primeSimilarityIdxPrimed(i_),1} = 'PrimedLowSimilarity';
                t.similarityCondition{end-5+primeSimilarityIdxControl(i_),1} = 'ControlLowSimilarity';
            elseif i_ == 3
                t.similarityCondition{end-10+primeSimilarityIdxPrimed(i_),1} = 'PrimedMidSimilarity';
                t.similarityCondition{end-5+primeSimilarityIdxControl(i_),1} = 'ControlMidSimilarity';
            end
        end
        

    end

end

tabTrials = struct2table(t);

writetable(tabTrials, 'trialDataPrimeTargetSimilarity.csv')






