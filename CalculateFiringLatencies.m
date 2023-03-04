clearvars -except sessions eeg datadir stimdir% do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

if ~exist('sessions', 'var')
    load([datadir filesep 'sessions.mat']);
end
tic
%% params
minResponsesPerUnit = 1; % do this for units firing to at least one image
segstartstop = [100 1000]; % for latency analyses
doplot = false; % dont plot when calling response_latency;

%% load responses
pr = load('priming_responses.bak.mat');

% get units that significantly respond to four or more stimuli in both
% conditions -> need to recalculate primign responses for this
respondingInBoth=(pr.consider_rs(:,:,1) & pr.consider_rs(:,:,2));
uidx =  sum(respondingInBoth, 2) >= minResponsesPerUnit;

%% loop over units with responses in both and calculate latency
nunits = sum(uidx);
aggregatelatencies = NaN(nunits,2); % median across effective stimuli
aggregatedurations = NaN(nunits,2); % median across effective
aggregatebslfr = NaN(nunits);

latlookup.region = {};
latlookup.isInterneuron = []; 
latlookup.bslfr = [];
latlookup.method = {};

fuidx = find(uidx);
disp(sprintf('doing latency analysis on %d units', nunits));

for u = 1:nunits
    
    fprintf('%d ', u);
    
    % get index to unit in sessions
    clusid = fuidx(u);    
    sessid = pr.cluster_lookup.sessid(clusid);
    channo = pr.cluster_lookup.channo(clusid);
    classno = pr.cluster_lookup.clusid(clusid);
    
    cherryno = find([sessions(sessid).cherries(:).channr] == channo & ...
        [sessions(sessid).cherries(:).classno] == classno);
    
    assert(numel(cherryno) == 1)
            
    % get stimuli indices
    ir = find(respondingInBoth(clusid, :));
    %ir = find(pr.consider_rs(clusid, :));
    stimnames = pr.stim_lookup(ir);
    
    % init output for p_burst
    BOB = [];
    EOB = [];
    SOB = [];
    
    % calculate baselinfiring across all trials in KHz
    bslfr = firing_rate(sessions(sessid).cherries(cherryno).trial, -500, 0)/1000;
    
    % loop over stimuli
    for si = 1:numel(stimnames)
        
        % get indices to trials of current stimulus
        tidx_stim = strcmp([stimnames{si} '.jpg'],...
                           sessions(sessid).condition.imagename);
        %loop over conditions
        for ci = 1:2
            
            % get indices of trials within 
            tidx_cond = sessions(sessid).condition.condition == ci;
            tidx = tidx_cond & tidx_stim;
            
            trialonsetts=sessions(sessid).condition.onset_time(find(tidx));
            trialspikets=sessions(sessid).cherries(cherryno).trial(tidx);
            allspikets=sessions(sessid).cherries(cherryno).allspiketimes;
            
            [latencyinfo BOB EOB SOB bslfr] = response_latency(trialspikets, ...
                                                  trialonsetts, ...
                                                  allspikets, ...
                                                  BOB, ...
                                                  EOB, ...
                                                  SOB, ...
                                                  segstartstop, ...
                                                  doplot, ...
                                                  bslfr);
            
            alllatencies{u,ci, si} = latencyinfo.medianlatency;
            alldurations{u,ci, si} = latencyinfo.medianduration;
            allbslfr{u,ci,si} = bslfr;
        end
    end
    
    for ci = 1:2
        aggregatelatencies(u,ci) = nanmedian([alllatencies{u,ci,:}]);
        aggregatedurations(u,ci) = nanmedian([alldurations{u,ci,:}]);
    end
    
    % save lookup vars
    latlookup.region{u} = pr.cluster_lookup.regionname{clusid};
    latlookup.bslfr(u) = nanmedian(nanmedian([allbslfr{u,:,:}]));
end

disp(sprintf('took %.2f hours', toc/3600));

save(sprintf('priming_latencies_min%d_resp.mat',minResponsesPerUnit), ...
     'aggregatelatencies', ...
     'allbslfr', ...
     'alllatencies', ...
     'alldurations', ...
     'aggregatedurations', ...
     'latlookup');

%% play chirp when done
load chirp.mat;
sound(y, Fs);