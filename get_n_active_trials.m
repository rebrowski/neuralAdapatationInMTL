function [n_active_trials ntrials]=get_n_active_trials(spikes, nspikes_for_active, ...
                                                      response_period)

%% function [n_active_trials ntrials]=get_n_active_trials(spikes, nspikes_for_active, response_period)
% is a helper function for e.g. binwise_signed_rank.m 
% takes cell or array of spiketimes (rows = trials) and returns the
% number of trials in which at least nspikes_for_active spikes are found

if ~exist('response_period') || isempty(response_period)
    response_period = [0 1000];
end

if ~exist('nspikes_for_active') || isempty(nspikes_for_active)
    nspikes_for_active = 1;
end

if ~iscell(spikes)
    spikes=grapes2cell(spikes);
end

n_active_trials = 0;
ntrials = length(spikes);

for t = 1: ntrials
    if sum(spikes{t} > response_period(1) & ...
           spikes{t} < response_period(2)) ...
            >= nspikes_for_active
        
        n_active_trials = n_active_trials + 1;
    end
end
