function [dispersion idx] = get_dispersion(latencies, n_values)
%% function [dispersion idx] = get_dispersion(latencies, n_values)
% computes a home-made version to quantify dispersion of spike
% onset latencies:
% takes a vector of onset-latencies and calculates the range of the
% n_values closest to the median
% returns NaN if the fewer than n_values are in the latencies
% vector
if ~exist('n_values', 'var')
    n_values = ceil(numel(latencies) * 2/3);
end

if ~isempty(latencies) && length(latencies)>= n_values
    
    % get absolut differences to the median of all values in latencies
    lat = abs(latencies - median(latencies));

    % sort values according absolute difference to median
    [sortmat, idx] = sortrows([lat; 1:length(lat); latencies]');
    
    % compute dispersion between the lowest and highest of the
    % n_values closest values to the median
    idx = sortmat(1:n_values,2);
    lat = latencies(idx);
    lat = sort(lat);
    dispersion = lat(end) - lat(1);

else
    dispersion = NaN;
    idx = NaN;
end
