clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

% this takes time as its ~7GB:
if ~exist('sessions', 'var')
    tic
    load([datadir, filsep, datafilename])
    toc/60
end

nsessions = numel(sessions);
%% which sessions to do 
sessionstodo = 1:nsessions

%% get the number of units across all sessions defined above
nclustotal = 0;
sessid = []; % lookup session for clusters
for s = sessionstodo
    nclustotal = nclustotal + numel(sessions(s).cherries);
    sessid = [sessid repmat(s, 1, numel(sessions(s).cherries))];
end

%% some basic params 
bsl_window = [-500 0]; %ms
response_window = [100 1000];%ms
preStimWindow = bsl_window; %ms -> To assess spreading activation

%% some basic infos about the design
nstim = 100; % 100 images were used
ncat = 10; % belonging to 10 categories
nstimpercat = nstim/ncat; % there were 10 images for each of the 10 categories

% do a lookup for stim/category
stim_lookup = cell(100,1);
cat_lookup = cell(100,1);
lc = 1;
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

%% init lookup and output
nconditions = 2; %1: primed, 2:control
zvals = NaN(nclustotal, nstim, nconditions); % zscores
zvalsPreStim = NaN(nclustotal, nstim, nconditions); % zscors for baseline primed vs. control

sessid = NaN(nclustotal,1); % index in session struct
subjid = NaN(nclustotal,1); % subject number (053, 023...)
channo = NaN(nclustotal,1); % channelnumber in session
clusid = NaN(nclustotal,1); % clusternumber on channel
clustype = cell(nclustotal,1); % mu or su
sitename = cell(nclustotal,1); % on which site was the
regionname = cell(nclustotal,1); % anatomical label as defined
hemisphere = cell(nclustotal,1); % lef or right


%% compute an Nunits X Nstimuli X Nconditions matrix with Z-values
cc = 1; % output counter
for s = sessionstodo
    nclus = numel(sessions(s).cherries);
    t0 = tic;
    disp(sprintf('#%d %s', s, sessions(s).name))
    
    for c=1:nclus
        
        [mfr sd] = firing_rate(sessions(s).cherries(c).trial, ...
                               bsl_window(1), bsl_window(2)); 
        % get z-score: note that when SD == 0 & frstim > 0, Z = Inf
        % lets take care of that somehow
        if sd == 0
            sd = 1;
        end
        
        for cat= 1:ncat
            for stim = 1:nstimpercat
                for cond = 1:2
                    % save indices of trials
                    idx = sessions(s).condition.stimulus == stim & ...
                          sessions(s).condition.category == cat & ...
                          sessions(s).condition.condition == cond;
                    id_ = (cat-1)*10 + stim;
                    
                    [frstim sdstim] = firing_rate(sessions(s).cherries(c).trial(idx), ...
                                         response_window(1), ...
                                         response_window(2));
                    
                    [frprestim sdprestim] = firing_rate(sessions(s).cherries(c).trial(idx), ...
                                         preStimWindow(1), ...
                                         preStimWindow(2));
                    
                    
                    zvals(cc, id_, cond) = (frstim - mfr)/sd;
                    zvalsPreStim(cc, id_, cond) = (frprestim - mfr)/sd;

                    if frstim == 0
                        frstim = 1;
                    end
                    
                    fr(cc, id_, cond) = frstim;
                    frPreStim(cc, id_, cond) = frprestim;
                end
            end
        end
        
        % save lookup for session, clus on chan , and site
        sessid(cc) =  s;
        subjid(cc) =  str2num(sessions(s).name(1:3));
        channo(cc) = sessions(s).cherries(c).channr;
        clusid(cc) = sessions(s).cherries(c).classno;
        sitename{cc} = sessions(s).cherries(c).site;
        clustype{cc} = sessions(s).cherries(c).kind;
        % lookup region, and hemisphere
        reg_ = 'other';
        for r = 1:numel(regions)              
            if sum(strcmp(sitename{cc}, regions(r).sites)) > 0
                reg_ = regions(r).name;
            end
        end
        regionname{cc} = reg_;
        hemisphere{cc} = 'L';
        if strfind(sitename{cc}, 'R') > 0
            hemisphere{cc} = 'R';
        end

        cc = cc + 1;
    end
    
end

%% save the resulting zvals and lookup
cluster_lookup = table(sessid, subjid, channo, clusid, sitename, ...
                       regionname, hemisphere, clustype);
save('zvals_condition.mat', '-v7.3', ...
     'cluster_lookup', 'stim_lookup', 'cat_lookup', 'regions', ...
     'zvals', 'zvalsPreStim', 'fr', 'frPreStim');
