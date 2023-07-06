function [latencyinfo BOB EOB SOB fr] = response_latency(trialspikets, ...
                                                  trialonsetts, ...
                                                  allspikets, ...
                                                  BOB, ...
                                                  EOB, ...
                                                  SOB, ...
                                                  segstartstop, ...
                                                  doplot, ...
                                                  bslfr, ...
                                                  bslfrmode)


%% function [latency burstonsets burstoffsets BOB EOB SOB] = response_latency(trialspikets, ...
%                                                  trialonsetts, ...
%                                                  allspikets, ...
%                                                  BOB, ...
%                                                  EOB, ...
%                                                  SOB, ...
%                                                  segstartstop, 
%                                                  doplot, ...
%                                                  bslfr, ...
%                                                  bslfrmode)
%
% computes the latency of firing corresponding to spikets in
% trialspikets. requires, all spikts of the cluster to compute
% poisson-burst detection (p_burst.m). To get burst on/offsets
% relative for each trial, trial-onset times are required. Returns
% a the median onset time of burst (latency) as well as cells with
% burst on/offsets in each trial, and the output of p_burst.
% assumes all timestamps in milliseconds.
% bslfrmode: either compute the baselinefiringrate in -500:0 ms
% relative stim onset (bslfrmode = 'bslIntervalOnly', default), or
% using all recorded spikes (bslfrmode = 'allspikes')
% note that if you choose 'bslIntervalOnly' and run the function
% over multiple conditions/stimuli per cell it is recommended to
% calculate baselinefiring across all trials outside of this
% function and provide this number as argument (bslfr)
% provide bslfr in Khz for it to work with p_burst later on



if ~iscell(trialspikets)
    trialspikets = grapes2cell(trialspikets);
end
ntrials = length(trialspikets);
assert(length(trialonsetts) == ntrials)
% consider burst in this period relative to trialonsetts
if ~exist('segstartstop', 'var') 
    segstartstop = [100, 1000]; % where to look for bursts
end

if ~exist('doplot', 'var')
    doplot = false;
end

% compute the average firing rate of the cell


if exist('bslfr', 'var') 
        
    fr = bslfr;
    %keyboard
else
        
    if ~exist('bslfrmode', 'var')
        bslfrmode = 'bslIntervalOnly';
    end

    bsl_fr = firing_rate(trialspikets, -500, 0)/1000; % baseline of trialspikets
    overall_fr = numel(allspikets)/(max(allspikets)-min(allspikets)); % allspikets
                                                                      % which of the above should be use?
                                                                      %fr = overall_fr;
    if strcmp(bslfrmode, 'bslIntervalOnly')
        fr = bsl_fr;
    else
        fr = overall_fr;
    end
end



method = '';

if fr*1000 > 2 
    method = 'poisson burst detection';
    %% higher than 2 Hz > use p_burst.m
    % if you have multiple groups of trials per channel, you might want
    % to do this only once per channel and pass BOB and his friends as arguments later
    % on as it takes a while
    if ~exist('BOB', 'var') || isempty(BOB)
        [BOB, EOB, SOB]=p_burst(allspikets, fr);
    end

    % get BOBs and EOB relative to trialonsetts for each trial
    BOBts = allspikets(BOB);
    EOBts = allspikets(EOB);

    for t = 1:ntrials
        burstonsets{t} = BOBts(BOBts > trialonsetts(t) + segstartstop(1) & ...
                               BOBts < trialonsetts(t) + segstartstop(2));
        
        % only check for offsets if there is at least one onset
        if ~isempty(burstonsets{t})
            firstonset = burstonsets{t}(1);
             burstoffsets{t} = EOBts(EOBts > firstonset & ...
                                EOBts < trialonsetts(t) + segstartstop(2));

        else
            burstoffsets{t} = [];
        end
        
        burstonsets{t} = burstonsets{t} - trialonsetts(t);
        burstoffsets{t} = burstoffsets{t} - trialonsetts(t);
    end
else 

    %% less than 2Hz firing rate -> get first spike after
    %%  segstartstop(1)
    method = sprintf('first spike %d ms after stimulus onset', segstartstop(1));
    for t = 1:ntrials
       burstonsets{t} = trialspikets{t}(find(trialspikets{t} > segstartstop(1) ...
                                        & trialspikets{t} < ...
                                             segstartstop(2), 1, ...
                                             'first'));
       burstoffsets{t} = burstonsets{t} + 1; % duration = 1ms
    end
    
    % they have to exist..
    BOB = [];
    EOB = [];
    SOB = [];
end

%% determine the latency, i.e., only get the first of the bursts in
%% a trial and calculate burst-duration where possible
latencies = [];

for t = 1:ntrials
    if ~isempty(burstonsets{t})
        latencies(t) = burstonsets{t}(1);
    else 
        latencies(t) = NaN;
    end    
end

%% get first offset
offsets = [];

for t = 1:ntrials
    if ~isempty(burstoffsets{t})
        offsets(t) = burstoffsets{t}(1);
    else 
        offsets(t) = NaN;
    end    
end

%% calculate durations form the above

durations = offsets - latencies;

nvalsclosesttomedian = ceil(sum(~isnan(latencies))*2/3); % tocompute dispersion

%% assign output
latencyinfo.latencies = latencies;
latencyinfo.medianlatency = nanmedian(latencies);
if ~isnan(latencyinfo.medianlatency)
    latencyinfo.dispersion = get_dispersion(latencies(~isnan(latencies)), ...
                                            nvalsclosesttomedian);
else
    latencyinfo.dispersion = NaN;
end

latencyinfo.nvalsclosesttomedian = nvalsclosesttomedian;
latencyinfo.method = method;
latencyinfo.burstonsets = burstonsets;
latencyinfo.burstoffsets = burstoffsets;
latencyinfo.burstdurations = durations;
latencyinfo.medianduration = nanmedian(durations);



%% plot things
if doplot
    figure
    x = -1000:2000;
    subplot(1,2,1);
    plot_raster(trialspikets,x,[],'k')
    title(sprintf('firing rate %.1f Hz, method: %s, latency:%.1f ms dispersion %.1f ms', ...
                  fr * 1000, method, latencyinfo.medianlatency, latencyinfo.dispersion));
    subplot(1,2,2);
    plot_raster(trialspikets,x,[],'k');
    plot_raster(burstonsets,x,[],'g');
    plot_raster(burstoffsets,x,[],'r');
    plot([segstartstop(1) segstartstop(1)], ylim, '-b')
    plot([segstartstop(2) segstartstop(2)], ylim, '-b')
end

