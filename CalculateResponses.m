clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

if ~exist('sessions', 'var')
    load([datadir, filesep, 'sessions.mat'])
end

load regions

nsessions = numel(sessions);
nregions = numel(regions);
%% which sessions to do 
sessionstodo = 1:nsessions

%% get the number of units across all sessions
nclustotal = 0;
sessid = []; % lookup session for clusters
for s = sessionstodo
    nclustotal = nclustotal + numel(sessions(s).cherries);
    sessid = [sessid repmat(s, 1, numel(sessions(s).cherries))];
end

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
%% get pvalues and lookuptables for analyses of interest
% paramters for binwise_signedrank and ranksum
bsl_period = [-500 0];
response_period = [0 1000];
binsize = 100;

pcrit = 0.001; % for the binwise_ranksum
pcrit_sr = .05;
pcrit_rs = 0.001;

min_fraction_of_active_trials = .5; % more than 5/10 trials have a spike within 0-1000 ms
consider_positive_only = true;
oldversion = false; % the original analysis has a inconsistency with
                   % plot_responses.m as it uses histc.m wrongly (see binwise_ranksum.m)
excludebelow2hz = true; % some version of the binwise ranksum test
                         % do not consider responses with fr < 2 Hz
                         % during the response period

% init output:
pvals_rs = NaN(nclustotal, nstim); % p < prcir & positive response
consider_rs = NaN(nclustotal,nstim);

% init lookup
subjid = NaN(nclustotal,1); % sub number
sessid = NaN(nclustotal,1); % index in session struct
channo = NaN(nclustotal,1); % channelnumber in session
clusid = NaN(nclustotal,1); % clusternumber on channel
sitename = cell(nclustotal,1); % on which site was the
regionname = cell(nclustotal,1); % anatomical label as defined
hemisphere = cell(nclustotal,1); % left or right
clustype = cell(nclustotal,1); % SU or MU

%% calculate pvals of responses per image
cc = 1;
for s = sessionstodo
    nclus = numel(sessions(s).cherries);
    t0 = tic;
    disp(sprintf('#%d %s', s, sessions(s).name))
    subj_ = str2num(sessions(s).name(1:3));
    for c = 1:nclus
        fprintf('%d ', c)
        distribution = [];
        for cat = 1:ncat
            for stim = 1:nstimpercat

                % save indices of trials
                idx = sessions(s).condition.stimulus == stim & ...
                      sessions(s).condition.category == cat & ...
                      sessions(s).condition.condition > 0;
                
                id_ = (cat-1)*10 + stim;

 
                % calculate binwise ranksum
                [p_ consider distribution] = ...
                    binwise_ranksum(sessions(s).cherries(c).trial(idx), ...
                                    sessions(s).cherries(c).trial, ...
                                    bsl_period, response_period, ...
                                    binsize,pcrit, min_fraction_of_active_trials, ...
                                    consider_positive_only, distribution, ...
                                    excludebelow2hz, oldversion);

                pvals_rs(cc,id_ ) = p_;
                consider_rs(cc,id_) = consider;
                

                clear consider p_
            end % stimulus
        end % category

        % save lookup for session, clus on chan , and site
        subjid(cc) = subj_;
        sessid(cc) = s;
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
    end % cluster   
    disp(sprintf(' - took %.1f minutes', toc(t0)/60));
end

%% plot n classes per session
for s = 1:nsessions
   nclasses(s) = numel(sessions(s).cherries);
end
fh = figure('color', 'w', 'Position', [200 200 800 300 ])
bar(nclasses)
xlabel('session #');
ylabel('N clusters');
legend(sprintf('N total = %d', sum(nclasses)));
legend boxoff
box off
print_png(['plots' filesep 'N_Clusters_per_Session.png'], fh);

%% get proportions of sign. responses per category and region and plot
nhemispheres = 2;
uhem = {'L', 'R'};
ncat = 10;
nstim = 100;
nclustotal = size(pvals_rs, 1);
nsig_rs = zeros(nhemispheres, nregions,ncat);
nsig_sr = zeros(nhemispheres, nregions,ncat);
ntot = zeros(nhemispheres, nregions,ncat);

tic
for cc = 1:nclustotal
    for stim = 1:nstim
        cat = floor((stim-1)/10)+1;
        for r = 1:nregions
            for h = 1:nhemispheres
                if strcmp(regions(r).name, regionname{cc}) && ...
                        strcmp(uhem{h}, hemisphere{cc})
                     if consider_rs(cc,stim) %pvals_rs(cc,stim) < pcrit_rs    
                        nsig_rs(h, r, cat) = nsig_rs(h,r,cat) + 1;                        
                    end
                    ntot(h,r,cat) = ntot(h,r,cat) + 1;
                end
            end % hemispheres
        end % regions
    end % category
end % clus
toc

%% save p values 
cluster_lookup = table(sessid, channo, clusid, sitename, regionname, ...
                       hemisphere, subjid, clustype);

save('category_responses.mat', '-v7.3', ...
     'cluster_lookup', ...
     'stim_lookup', ...
     'cat_lookup', ...
     'regions', ...
     'pvals_rs', ...
     'consider_rs', ...
     'pcrit_rs', ...
     'nsig_rs', ...
     'ntot');


