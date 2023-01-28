clearvars -except sessions eeg

%% load data
eegname = '~/projects/ospr/secondlevel/ieeg_linked_mastoids_256Hz.mat'; % adjust this to your situation

if ~exist('eeg', 'var')
    load(eegname); %loads var eeg into the workspace
end
nsessions = numel(eeg);

%% go through ERPS and get latencies of first and second peak for primed and control condition
%% time-windows
% regions & time window parameters
regions(1).ieegsites = {'AL', 'AR'};
regions(1).sites =     {'LA', 'RA'};
regions(1).name = 'AM';
regions(1).tw1from = 200; % ms 
regions(1).tw1to = 400;
regions(1).tw2from = 400;
regions(1).tw2to = 750;

regions(2).ieegsites = {'AHL', 'AHR', 'MHL', 'MHR', 'PHL', 'PHR'};
regions(2).sites =     {'LAH', 'RAH', 'LMH', 'RMH', 'LPH', 'RPH'};
regions(2).name = 'HC';
regions(2).tw1from = 200; % ms 
regions(2).tw1to = 400;
regions(2).tw2from = 400;
regions(2).tw2to = 750;

regions(3).ieegsites = {'ECL', 'ECR'};
regions(3).sites =     {'LEC', 'REC'};
regions(3).name = 'EC';
regions(3).tw1from = 200; % ms 
regions(3).tw1to = 400;
regions(3).tw2from = 400;
regions(3).tw2to = 750;

regions(4).ieegsites = {'PHCL', 'PHCR'};
regions(4).sites =     {'LPHC', 'RPHC'};
regions(4).name = 'PHC';
regions(4).tw1from = 200; % ms 
regions(4).tw1to = 400;
regions(4).tw2from = 400;
regions(4).tw2to = 750;

% aggregate data 
for sessi = 1:nsessions

    for regi = 1:numel(regions)
          % get indices to samples for timewindows
        tw1idx = eeg(sessi).stime > regions(regi).tw1from & eeg(sessi).stime <= regions(regi).tw1to;
        tw2idx = eeg(sessi).stime > regions(regi).tw2from & eeg(sessi).stime <= regions(regi).tw2to;

        % 1. get channels matching region
        nsites_in_region = numel(regions(regi).ieegsites);
        siteidx = false(1,eeg(sessi).nsites);
        for sitei = 1:nsites_in_region
            i_ = find(strcmp(regions(regi).sites(sitei), ...
                eeg(sessi).sites));
            if ~isempty(i_)
                siteidx(i_) = true;
            end

        end

        if sum(siteidx) > 0
            siteidxf = find(siteidx);
            clear dat

            for chani = 1:numel(siteidxf)

                for cond = 1:2

                    tidx = eeg(sessi).condition.condition == cond & ...
                        ~eeg(sessi).isartefact(siteidxf(chani),:);

                    % 1. average over trials per channel
                    dat(chani, cond, :) ...
                        = squeeze(mean(eeg(sessi).sdata(siteidxf(chani),tidx,:),2));

                end % conditions
            end % channels in site

            % 2. average over channels, get peak-latencies in timewindow

            for cond = 1:2
                % erp
                erp(regi).dat(sessi,cond, :) = ...
                    squeeze(mean(dat(:,cond,:),1) .* -1); % *. -1 as in the figure to reflect actual polarity

                
                % get maxima and minima and latencies in TW1
                [m_ i_] = max(squeeze(erp(regi).dat(sessi, cond,tw1idx)));
                tax_ = eeg(regi).stime(tw1idx);
                erp(regi).tw1max(sessi, cond) = m_ ;
                erp(regi).tw1maxlatency(sessi, cond) = tax_(i_);
               

                %erp in timewindow> 
                %plot(tax_, squeeze(erp(regi).dat(sessi, cond,tw1idx)))
                %hold on
                %plot(tax_(i_), m_, 'or');
               
                [m_ i_] = min(squeeze(erp(regi).dat(sessi, cond,tw1idx)));
                erp(regi).tw1min(sessi, cond) =  m_;
                erp(regi).tw1minlatency(sessi,cond) = tax_(i_);
              
                % get maxima and minima and latencies in TW2
                [m_ i_] = max(squeeze(erp(regi).dat(sessi, cond,tw2idx)));
                tax_ = eeg(regi).stime(tw2idx);
                erp(regi).tw2max(sessi, cond) = m_ ;
                erp(regi).tw2maxlatency(sessi, cond) = tax_(i_);
               
               
                [m_ i_] = min(squeeze(erp(regi).dat(sessi, cond,tw2idx)));
                erp(regi).tw2min(sessi, cond) =  m_;
                erp(regi).tw2minlatency(sessi,cond) = tax_(i_);
            end

            erp(regi).subject(sessi) = eeg(sessi).condition.subject;
        end % if


    end % regions
    erp(regi).regionname = regions(regi).name;

    
end % sessions

%% print results
for r = 1:numel(regions)

    sessidx = erp(r).tw2max(:,1) ~= 0;
    disp(' ');
    disp(sprintf('%s, N=%d', regions(r).name, sum(sessidx)));
    disp('=============')
%     
    % TW1 minimum/ negative peak
    signrankP = signrank(erp(r).tw1minlatency(sessidx, 1), erp(1).tw1minlatency(sessidx, 2));
    [h ttestp] = ttest(erp(r).tw1minlatency(sessidx, 1), erp(1).tw1minlatency(sessidx, 2));
    disp(sprintf('Min %d-%dms: med(iqr)_{primed}=%.2f(%.2f), med(iqr)_{control}=%.2f(%.2f), p_{signrank}=%.3g, p_{ttest}=%.3g', ...
        regions(r).tw1from, regions(r).tw1to, ...
        median(erp(r).tw1minlatency(sessidx,1)),iqr(erp(r).tw1minlatency(sessidx,1)), ...
        median(erp(r).tw1minlatency(sessidx,2)),iqr(erp(r).tw1minlatency(sessidx,2)), ...
        signrankP, ttestp));    
    
    %TW2, maximum/positive peak
    signrankP = signrank(erp(r).tw2maxlatency(sessidx, 1), erp(1).tw2maxlatency(sessidx, 2));
    [h ttestp] = ttest(erp(r).tw2maxlatency(sessidx, 1), erp(1).tw2maxlatency(sessidx, 2));
    disp(sprintf('Max %d-%dms: med(iqr)_{primed}=%.2f(%.2f), med(iqr)_{control}=%.2f(%.2f), p_{signrank}=%.3g, p_{ttest}=%.3g', ...
        regions(r).tw2from, regions(r).tw2to, ...
        median(erp(r).tw2maxlatency(sessidx,1)),iqr(erp(r).tw2maxlatency(sessidx,1)), ...
        median(erp(r).tw2maxlatency(sessidx,2)),iqr(erp(r).tw2maxlatency(sessidx,2)), ...
        signrankP, ttestp));
    
end


