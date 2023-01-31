%% this scripts calculates 'tuning curves' for each unit with at
%% least N response-eliciting stimuli

clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

%% definition of recording sites in anatomical regions
regions(1).sites =     {'LA', 'RA'};
regions(1).name = 'AM';

regions(2).sites =     {'LAH', 'RAH', 'LMH', 'RMH', 'LPH', 'RPH', 'LEC', 'REC', 'LPHC', 'RPHC'};
regions(2).name = 'HC, EC, PHC';
nregions = numel(regions);

%% load response-criterion
cr = load('category_responses.mat');

%% load zscores and firing rates during primed vs. control condition
fr = load('zvals_condition.mat');

%% do this for units responding to at least 4 stimuli
min_n_re_stim = 4

    %% get indices of units with at least min_n_re_stims
    uidx = sum(cr.consider_rs,2) >= min_n_re_stim; % p/c combined;
    nunits = sum(uidx);
    cluster_infos = fr.cluster_lookup(uidx,:);
    nsu = sum(strcmp('SU', cluster_infos{:,'clustype'}) & ~strcmp('other', cluster_infos{:,'regionname'}));
    uidx = find(uidx); % get the actual indices
    
    % init output
    nconditions = 2; % primed = 1, control = 2;
    nstimuli = 100;
    tc = NaN(nunits, 100, nconditions); % normalized to max of control condition
    tcfr = NaN(nunits, 100, nconditions); % non-normalized fr
    
    % control condition
    tcfr = NaN(nunits, 100, nconditions); % raw-firing rates
    isresponse = NaN(nunits, 100); % save wheter it is response or not
    
    for u = 1:nunits
        % 1. sort according to mean firing rate of primed/control
        [rmagnitude sidx] = sort(mean(fr.fr(uidx(u),:,:),3), ...
            'descend');
        tcfr(u,1:100, :) = fr.fr(uidx(u), sidx, :);
        tc(u,1:100,:) = fr.fr(uidx(u), sidx, :);

        % 3. apply sorting to primed, control condition and response-criterion
        is_response(u, :) = cr.consider_rs(uidx(u),sidx);

        % 4. normalize
        max_ = tc(u,1,2); % max fr control
        tc(u, 1:100, :) = tc(u, 1:100, :)./max_;
    end

    %% save the tuning curves for further stuff
    fname= sprintf('tuningCurvesMin%dResponsesPerUnit', ...
                                min_n_re_stim);
    save(fname, 'tc', 'tcfr', 'is_response', 'cluster_infos');