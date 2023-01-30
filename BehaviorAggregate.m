
clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

%% this script loads the raw data (sessions.mat) and aggregates relevant behavioral data 
if ~exist('sessions', 'var')
    load([datadir filesep 'sessions.mat']); % this may take a while as sessions.mat is 4.7 GB
end

%% Basic Parameters
manmadebutton = 1; % arrow key to the left
naturalbutton = 2; % arrow key to the right
ncat = 10; % number of categories
ucond = [1 2]; % category repeated = 1, non-repeated = 2;
ncondition = 2;
nstim = 100; % 100 images were used
nstimpercat = nstim/ncat; % there were 10 images for each of the 10 categories
rtout = 2.5; % Zval threshold for outlier removal

%% do a lookup for stim/category
stim_lookup = cell(100,1);
cat_lookup = cell(100,1);
lc = 1; % counter for long data
for cat = 1:ncat
    for stim = 1:nstimpercat
        id = find(sessions(1).condition.stimulus == stim & ...
                  sessions(1).condition.category == cat, 1, 'first');
        imname = sessions(1).condition.imagename{id};
        
        stim_lookup{lc} = imname(1:strfind(imname, '.jpg')-1);
        cat_lookup{lc} = imname(1:regexp(imname, '_[1-9]')-1);
        lc = lc + 1;
    end
end
stimlookup = table(stim_lookup, cat_lookup);
categorynames = cat_lookup(1:10:end);
conditionnames = {'primed', 'control'};

%% Initialize some variables for aggregate values
% percentage correct
pcorrect = NaN(1,length(sessions)); % overall
pcorrcat_wide = NaN(length(sessions), ncat, ncondition);
rtcat_wide = NaN(length(sessions), ncat, ncondition);
subid_wide = NaN(1,length(sessions));

pcorrcat_long = NaN(1, length(sessions) * ncat * ncondition);
rtcat_long = NaN(1, length(sessions) * ncat * ncondition);
subid_long = NaN(1, length(sessions) * ncat * ncondition);
condition_long = NaN(1,length(sessions) * ncat * ncondition);
conditionname_long = cell(1, length(sessions) * ncat * ncondition);
categoryname_long = cell(1, length(sessions) * ncat * ncondition);


% restrict to trials where in control condition the previous
% stimulus was form the same meta-category (manmade or natural),
% i.e. response priming
pcorrcat_wide_sr = NaN(length(sessions), ncat, ncondition);
rtcat_wide_sr = NaN(length(sessions), ncat, ncondition);

pcorrcat_long_sr = NaN(1, length(sessions) * ncat * ncondition);
rtcat_long_sr = NaN(1, length(sessions) * ncat * ncondition);
% $$$ subid_long_sr = NaN(1, length(sessions) * ncat * ncondition);
% $$$ condition_long_sr = NaN(1,length(sessions) * ncat * ncondition );
% $$$ conditionname_long_sr = cell(1, length(sessions) * ncat * ncondition);
% $$$ categoryname_long_sr = cell(1, length(sessions) * ncat * ncondition);

%% Aggregate data

lc = 1;

for s = 1:length(sessions)
    % subid 
    subid = str2num(sessions(s).name(1:3));
    
    % get indices of manmade and natrual things
    mmidx = sessions(s).condition.category > 5;

    nidx = sessions(s).condition.category < 6;
    ntrials = sum(mmidx) + sum(nidx);
    assert(ntrials == 1010);
    
    for t=1:ntrials
        
        
        % code correct/incorrect reponses
        if sessions(s).condition.category(t) <= 5 
            sessions(s).condition.correct_response(t) = 2; % natural things (right button)
        else
            sessions(s).condition.correct_response(t) = 1; % manmade things (left button)
        end
        
        if sessions(s).condition.response(t) == sessions(s).condition.correct_response(t)
            sessions(s).condition.correct(t) = 1;
        else
            sessions(s).condition.correct(t) = 0;
        end
        
        % code supercategory repeitions
        if t > 1 && sessions(s).condition.condition(t-1) ~= 0 % not the first trial in a run 
            if sessions(s).condition.correct_response(t) == ...
                    sessions(s).condition.correct_response(t-1)
                sessions(s).condition.super_category_repetition(t) = 1; % for primed condition, this is always the case
            else
                sessions(s).condition.super_category_repetition(t) = 2; 
            end
        else
            sessions(s).condition.super_category_repetition(t) = 0;
        end
    end
    
    % store accuracy per session
    pcorrect(s) = ...
        (sum(sessions(s).condition.response(mmidx) == manmadebutton) + ...
         sum(sessions(s).condition.response(nidx) == naturalbutton))/...
         ntrials;
    
    pcmm(s) = sum(sessions(s).condition.response(mmidx) == manmadebutton)./sum(mmidx);
    pcn(s) = sum(sessions(s).condition.response(nidx) == naturalbutton)./sum(nidx);
    
    % get median rts for correct
    zrt = zscore(sessions(s).condition.rt); 
    rtinidx = zrt > -1*rtout & zrt < rtout;

    rtoutlierspercent(s) = sum(rtinidx)/ntrials;
    rtmm(s) = mean(sessions(s).condition.rt(mmidx & rtinidx));
    rtn(s) = mean(sessions(s).condition.rt(nidx & rtinidx));
    rtoverall(s) = mean(sessions(s).condition.rt((nidx|mmidx) & rtinidx));
    
    % accuracy&RT per condition
    for cond = 1:2 % primed vs. control
                        
            % indices current condition
            condidx = sessions(s).condition.condition == cond;
            
             pcorrect_condition(s, cond) = ...
        (sum(sessions(s).condition.response(mmidx & condidx) == manmadebutton) + ...
         sum(sessions(s).condition.response(nidx & condidx) == naturalbutton))/...
         ntrials*2;
            
            rtcorrect_condition(s, cond) = ...
                mean(sessions(s).condition.rt(condidx & rtinidx));
    end
      
    
    % do rt/%correct per category & condition
    subid_wide(s) = subid;

    for c = 1:ncat
        
        % indices to current category
        cidx = sessions(s).condition.category == c;
        % indices to correct responses
        corrbutton = naturalbutton;
        if c > 5
            corrbutton = manmadebutton;
        end
        corridx =sessions(s).condition.response ...
                                 == corrbutton; 


        for cond = 1:2 % primed vs. control
                        
            % indices current condition
            condidx = sessions(s).condition.condition == cond;
                       
            % correct responses to current category/condition combi
            ccidx = condidx & cidx & corridx;  
                                    
            pc = sum(ccidx)./ sum(cidx & condidx); % 10 stim/categorz 5 trials/condition
            rt = sessions(s).condition.rt(ccidx & rtinidx);

            % wide
            pcorrcat_wide(s,c,cond) = pc; 
            rtcat_wide(s,c,cond) = mean(rt);
            % long
            pcorrcat_long(lc) = pc;
            rtcat_long(lc) = mean(rt);
            subid_long(lc) = subid;
            catid_long(lc) = c;
            categoryname_long{lc} = categorynames{c};
            condid_long(lc) = cond;
            conditionname_long{lc} = conditionnames{cond};            
            
            % only consier trials with meta-category-repeition in
            % control condition (i.e, remove response priming
            % effect)

            % TR 28.7.2021 - control targets primed by same metacategory
            % (this is true for primed condition anyway)
            withinSameResponseIdx = ...
                sessions(s).condition.super_category_repetition == 1;
                        
            % get the appropriate indices
            ccidxmr = ccidx & withinSameResponseIdx;
            pc = sum(ccidxmr)./ sum(cidx & condidx & withinSameResponseIdx); % 10 stim/categorz 5 trials/condition
            rt = sessions(s).condition.rt(ccidx & rtinidx & withinSameResponseIdx);
            
            % wide
            pcorrcat_wide_sr(s,c,cond) = pc; 
            rtcat_wide_sr(s,c,cond) = mean(rt);
            % long
            pcorrcat_long_sr(lc) = pc;
            rtcat_long_sr(lc) = mean(rt);

            lc = lc + 1;
            
        end
    end
end

% save some aggregated data to disk
save('reactiontimes_primed_control_category.mat', 'rtcat_wide_sr', 'pcorrcat_wide_sr', 'subid_wide', 'rtcat_wide', 'pcorrcat_wide');



