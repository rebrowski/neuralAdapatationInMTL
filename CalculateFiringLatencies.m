%% this scripts calculates for each neuronal response according the
% binwise ranksum test a latency value for pthe primed and control
% condition

clearvars -except sessions

if ~exist('sessions', 'var')
    load sessions;
end
tic
%% params
minResponsesPerUnit = 4; % in both primed & control condition
segstartstop = [100 1000]; % for latency analyses
doplot = false; % dont plot when calling response_latency;

%% load binwise ranksum results per condition:
%pr = load('priming_responses.mat');

%% load binwise signed rank for p/c collapsed
cr = load('category_responses.mat');

%% load spikewidht/ pyrInt classification 
sw = load('ospr_spikewidth.mat');

% get units that significantly respond to primed AND control trials
% respondingInBoth=(pr.consider_rs(:,:,1) & pr.consider_rs(:,:,2));
% uidx =  sum(respondingInBoth, 2) >= minResponsesPerUnit;

uidx = 

disp(sprintf(['found %d repsonses that in units least %d responses ' ...
              'which are significant during both, primed & control condition'], sum(sum(respondingInBoth)), minResponsesPerUnit));
disp(sprintf(['these responses are found in  %d units, i.e., %.2f responses per unit on average '], ...
             sum(uidx), sum(sum(respondingInBoth))/sum(uidx)));

%% loop over units with responses in both and calculate latency
nunits = sum(uidx);
aggregatelatencies = NaN(nunits,2); % median across effective stimuli
aggregatedurations = NaN(nunits,2); % median across effective
aggregatebslfr = NaN(nunits);

latlookup.region = {};
latlookup.isInterneuron = []; 

% celltype_ison_crit2 = (spikewidth > 0.6) + 1;
% 1: interneurons
% 2: principal cells, 
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
    stimnames = pr.stim_lookup(ir);
    
    % init BOB and his friends 
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
    latlookup.isInterneuron(u) = sw.celltype_ison_crit2(clusid);
    latlookup.bslfr(u) = nanmedian(nanmedian([allbslfr{u,:,:}]));
end

disp(sprintf('took %.2f Minutes', toc/60));

save(sprintf('priming_latencies_min%dresponsesperunit.mat', minResponsesPerUnit), ...
     'aggregatelatencies', ...
     'allbslfr', ...
     'alllatencies', ...
     'alldurations', ...
     'aggregatedurations', ...
     'latlookup');

%% play chirp when done
load chirp.mat;
sound(y, Fs);