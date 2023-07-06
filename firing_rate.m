function [fr sd count minc maxc]= firing_rate(trial, from, to)

%% function [fr sd count] = firing_rate(trial, from, to)
% compute the mean (sd: standard deviation) firing rate (fr) in Hz for a cell array with spike-times 
% trial{1} = [-301, -277, ...], trial{2} = [-113 2 22 565...]
% from and to are times in ms denoting a period during which firing
% rate should be returned (e.g. for baseline-firing from = -400, to
% = 0 (ms)

if ~iscell(trial)
    trial = grapes2cell(trial);
end

duration = abs(to - from) / 1000; % seconds

if ~iscell(trial)
    spikes=grapes2cell(trial);
end

ntrials = length(trial);

for t = 1:ntrials
    count(t) = sum(trial{t} > from & ...
                    trial{t} < to);
end

minc = min(count);
maxc = max(count);
sd = std(count./duration);
fr = sum(count)/ duration / ntrials;

