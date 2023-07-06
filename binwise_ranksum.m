function [pval, consider, distribution] = binwise_ranksum(spikes, bsl_spikes, bsl_period, response_period, ...
                                               binsize, pcrit, ...
                                               min_fraction_of_active_trials, ...
                                               consider_positive_only, ...
                                                      distribution, ...
                                                      excludebelow2hz, oldversion)
%% function [pval, consider, distribution] = binwise_ranksum(spikes, bsl_spikes, bsl_period, response_period, ...
%                                               binsize, pcrit, ...
%                                               min_fraction_of_active_trials, ...
%                                               consider_positive_only, ...
%                                                   distribution, ...
%                                                      excludebelow2hz, oldversion))
%
% A.K.A the screening criterion

if ~exist('response_period') || isempty(response_period)
    response_period = [0 1000]; % milliseconds
end
if ~exist('bsl_period') || isempty(bsl_period)
    bsl_period = [-500 0]; % milliseconds
end
if ~exist('binsize', 'var') || isempty(binsize)
    binsize = 100; % milliseconds
end
if ~exist('pcrit', 'var') || isempty(pcrit)
    pcrit = 0.05;
end
if ~exist('min_fraction_of_active_trials', 'var') || isempty(min_fraction_of_active_trials)
    min_fraction_of_active_trials = 0.5;
end
if ~exist('consider_positive_only', 'var') || isempty(consider_positive_only)
    consider_positive_only = true;
end
if ~exist('excludebelow2hz', 'var') || isempty(excludebelow2hz)
    % original version in plot_responses.m does not consider
    % responses if firing rate during response-period is below 2 Hz
    excludebelow2hz = false;
end

if ~exist('oldversion', 'var') || isempty(oldversion)
    oldversion = false;
    % before Apr. 13 2018, we neglected the fact that histc.m
    % returns as many values as there edges -the last value
    % reflects to count of values matching the last value in Edges
    % this funtion did include it to compute the mean across bins
    % in the baseline-distrubtion - to reproduce this behavior, set
    % oldversion to true;
end

tmin = response_period(1); %ms
tmax = response_period(2);
tmin_baseline = bsl_period(1);
tmax_baseline = bsl_period(2);

baseline_duration = abs(tmin_baseline - tmax_baseline) / 1000; % seconds
response_period_duration = binsize / 1000; % seconds

if mod(baseline_duration * 1000, binsize) ~= 0 || ...
     mod(response_period_duration * 1000, binsize) ~= 0
     error(['Choose values for response- and baseline period which are ' ...
            'integer multiples of binsize'])
end

if ~iscell(spikes)
    spikes=grapes2cell(spikes);
end

if ~iscell(bsl_spikes)
    bsl_spikes=grapes2cell(bsl_spikes);
end

% get bsl distribution of spike counts
if ~exist('distribution', 'var') || isempty(distribution)

    if oldversion
        %% old and slightly wrong
        distribution = cellfun(@(x) mean(histc(x,tmin_baseline:binsize: ...
                                               tmax_baseline)), bsl_spikes, ...
                               'UniformOutput', false);
        
        for i = 1:numel(distribution)
            if isempty(distribution{i})
                distribution{i} = 0;
            end
        end
        distribution = [distribution{:}];

    else
        %% new and consistent with plot_responses.m
        distribution = zeros(1, numel(bsl_spikes));
        for b = 1:numel(bsl_spikes)
            h_ = histc(bsl_spikes{b},...
                       tmin_baseline:binsize: tmax_baseline);
            
            % 
            distribution(b) = mean(h_(1:end-1));
            
        end
    end
end

ntrials = length(spikes);

for t = 1:ntrials
    
    % generate component histograms for response period
    N = histc(spikes{t}, tmin:binsize:tmax);
    if isempty(N)
        N = zeros(1, numel(tmin:binsize:tmax));
    end
    
    N(end-1) = N(end-1) + N(end); % last bin is the count of x == tmax
                                  % (matlab R2012b...)
    N = N(1:end-1); 
    
    N2 = histc(spikes{t},tmin+binsize/2:binsize:tmax);
    if isempty(N2)
        N2 = zeros(1, numel(tmin+binsize/2:binsize:tmax));
    end
    
    N2(end-1) = N2(end-1) + N2(end);
    N2 = N2(1:end-1);
    
    % combine overapping histograms into one
    n_ = length(N) + length(N2);
    M(1:2:n_) = N;
    M(2:2:n_) = N2;
    
    % normalize counts not necessary anymore as the mean across binned
    % baseline is taken above
    % M = M ./ response_period_duration;
    
    % save
    rs_bins(t,:) = M;
end

% do ranksum per bin
nbins = length(rs_bins(1,:));
distribution_median = median(distribution);
for b = 1:nbins
    
    if oldversion
        direction(b) = sign(mean(rs_bins(:, b)) - distribution_median);
    else
        %% according to florians plot_responses.m:
        direction(b) = sign(median(rs_bins(:, b)) - distribution_median);
    end
    [p, h] = ranksum(rs_bins(:, b), distribution);
    pvals(b) = p; 
end

% correction for multiple comparisons (simes procedure)
[pvals, idx] = sort(pvals);
direction = direction(idx);
pvals_sime = nbins.*pvals./[1:nbins];
[pvals_sime, idx] = sort(pvals_sime);
pvals_sime(pvals_sime > 1) = 1;
direction = direction(idx);

% is p < crit & positive?
if consider_positive_only
    consider = sum(direction == 1 & pvals_sime < pcrit) > 0;
    p_index = find(pvals_sime < pcrit & direction ...
                   == 1, 1,'first');
    if ~isempty(p_index)
        pval = pvals_sime(p_index);
    else % not significant and positive
        p_index = find(direction == 1, 1, 'first');
        if ~isempty(p_index) % report lowest p-val of a
                             % postive bin
            pval = pvals_sime(p_index);
        else % no positive bins
            pval = 1;
        end
    end
else
    consider = sum(pvals_sime < pcrit) > 0;
    pval = pvals_sime(1);
end

n_active_trials = get_n_active_trials(spikes);
consider = consider * (n_active_trials/ntrials > ...
                       min_fraction_of_active_trials);


fr = firing_rate(spikes, response_period(1), response_period(2));
if fr < 2 & excludebelow2hz
    consider = false;
end


