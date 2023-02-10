 %% this script loads segmented EEG data an plots Figure 1C

clearvars -except sessions eeg datadir % do not clear the large variables if they are already loaded as this takes a few minutes
%% load things, set paths if necessary
if ~exist('datadir', 'var')
    startup
end

%% path to data % load
eegname = [datadir filesep 'ieeg_linked_mastoids_256Hz.mat']; % adjust this to your situation
if ~exist('eeg', 'var')
    load(eegname); %loads var eeg into the workspace
end

%% some params for statistical procedures, filenames, etc.
nsessions = numel(eeg);
alpha = 0.01;
clusteralpha = 0.001;
nperms = 1000;


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
figh.PaperPosition = [0 0 4.2 9.2];
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
    
    % TBD: store data in table for ncomms    
    timeax = eeg(1).stime(fromi:toi);
    datprimed = squeeze(erp(regi).dat(:,1,fromi:toi)) .* -1;
    datcontrol = squeeze(erp(regi).dat(:,2,fromi:toi)) .* -1;
 
    %condition = [repmat({'primed ERP'}, size(datprimed,1), 1); repmat({'control ERP'}, size(datcontrol,1),1)]   
    t = array2table([datprimed; datcontrol]);
    t.condition = [repmat({'primed ERP'}, size(datprimed,1), 1); repmat({'control ERP'}, size(datcontrol,1),1)];
    vl = strcat(cellstr(num2str(floor(timeax)'))', ' ms');
    t.Properties.VariableNames = [vl, 'condition'];
    tablename = ['source_data_files_ncomms', filesep, 'SourceDataFigure1C.xlsx'];
    writetable(t, tablename, 'Sheet', regions(regi).name); 
    
    if regi > 1
        set(handles(end), 'visible','off'); % this is the legend
    else
        set(handles(end), 'Location','best')
    end
    colormap(cbhandle, flipud(hot));
    cbh = colorbar(cbhandle);
    cbh.Position(1) = cbh.Position(1) + 0.08;
    cbh.Ticks = [0 0.2];
    cbh.TickLabels = {'0', '0.2'};
    cbh.Label.String = 'p';
    cbh.Label.Rotation = 0;
    cbh.Label.Position = [0.5 0.5 0];
    lbh = xlabel(handles(1), 'milliseconds');
    %if regi == 4 
    lbh.Position(2)= lbh.Position(2)-6;
    ylabel(handles(1), 'ERP amplitude (mV)');
    title(handles(1),regions(regi).name);
end

set(findall(gcf,'type','text'),'FontSize',fontSize);%,'fontWeight','normal');
set(findall(0,'type','axes'),'FontSize',fontSize)%,'fontWeight','normal');
print(figh,['plots', filesep, 'Figure1C'], '-dpng', '-r600')

