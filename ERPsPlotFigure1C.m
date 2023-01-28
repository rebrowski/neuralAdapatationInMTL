%% this script loads segmented EEG data an plots Figure 1C

clearvars -except sessions eeg % do not clear the large variables if they are already loaded as this takes a few minutes

%% path to data % load
eegname = '~/projects/ospr/secondlevel/ieeg_linked_mastoids_256Hz.mat'; % adjust this to your situation
if ~exist('eeg', 'var')
    load(eegname); %loads var eeg into the workspace
end

%% some params for statistical procedures, filenames, etc.
nsessions = numel(eeg);
alpha = 0.01;
clusteralpha = 0.001;
nperms = 1000;
erpoutfn = 'ERPsFigure';

%% definition of recording sites in anatomical regions
regions(1).ieegsites = {'AL', 'AR'};% Amygdalae
regions(1).sites =     {'LA', 'RA'};
regions(1).name = 'AM';

regions(2).ieegsites = {'AHL', 'AHR', 'MHL', 'MHR', 'PHL', 'PHR'}; % Hippocampi
regions(2).sites =     {'LAH', 'RAH', 'LMH', 'RMH', 'LPH', 'RPH'};
regions(2).name = 'HC';

regions(3).ieegsites = {'ECL', 'ECR'}; % entorhinal cortices
regions(3).sites =     {'LEC', 'REC'};
regions(3).name = 'EC';

regions(4).ieegsites = {'PHCL', 'PHCR'}; % parahippocampal cortices
regions(4).sites =     {'LPHC', 'RPHC'};
regions(4).name = 'PHC';

nregions = numel(regions);

%% aggregate ERPs by anatomical regions
for regi = 1:nregions
    
    sesscount = 1;
    
    for sessi = 1:nsessions
        
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
            end % cahnnels in site
            
            % 2. average over channels 
            for cond = 1:2
               erp(regi).dat(sesscount,cond, :) = ...
                   squeeze(mean(dat(:,cond,:),1));
            end
            sesscount = sesscount + 1;
        end % if 
        
        
    end % sessions
end % regions

%% Do Plot and Stats 
condcolors = [0 0 1; 1 0 0];
from = -500; % ms
to = 1500; % ms

fromi = find(eeg(1).stime < from, 1, 'last');
toi = find(eeg(1).stime > to, 1, 'first' );

figh = figure('color', 'w', 'visible', 'on');
figh.PaperUnits = 'inches';
figh.PaperPosition = [0 0 3.7 9.2];
% display it somewhat similar to what will be plotted
figh.Position = [200 200  figh.PaperPosition(3)*150 figh.PaperPosition(4)*150];
fontSize = 8;

for regi = 1:nregions
    subplot(4,1,regi)

    [handles k_ cbhandle] = ...
        perm_multcompare_pairwiseOnly(squeeze(erp(regi).dat(:,1,fromi:toi)) .* -1,...
                         squeeze(erp(regi).dat(:,2,fromi:toi)) .* -1, ...
               nperms,clusteralpha, alpha, eeg(1).stime(fromi:toi), condcolors, {'primed', ...
                        'control'}, gca, 5);

    if regi > 1
        set(handles(end), 'visible','off'); % this is the legend
    else
        set(handles(end), 'Location','best')
    end
    colormap(cbhandle, flipud(hot));
    xlabel(handles(1), '');
    ylabel(handles(1), 'Amplitude (mV)');
    title(handles(1),regions(regi).name);
end

set(findall(gcf,'type','text'),'FontSize',fontSize);%,'fontWeight','normal');
set(findall(0,'type','axes'),'FontSize',fontSize)%,'fontWeight','normal');

print(figh, erpoutfn, '-dpng', '-r600')

