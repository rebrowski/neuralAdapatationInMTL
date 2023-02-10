%% plot figures like 2A 

clearvars -except sessions eeg datadir stimdir

%% params 
cellids = [3628; 900; 4236; 3943; 900]; % for what units should the figures be done?
howmany = 16; %stimuli per cell

if ~exist('stimdir', 'var')
    startup
end
if ~exist('sessions', 'var')
    load([datadir filesep 'sessions.mat']);
end

%% load lookups
cr = load('category_responses.mat');
zc = load('zvals_condition');

%% loop through cells to plot
for cid = 1%:numel(cellids)
    cellid = cellids(cid);

    outfn = sprintf('%s%sospr_su_example_%d_%s', 'plots', filesep, cellid, cr.cluster_lookup{cellid,'sitename'}{1}); 

    figh = figure('color', 'w', 'visible', 'off');
    figh.PaperUnits = 'inches';
    figh.PaperPosition = [0 0 7.4 6];
          
    %% get info on unit to plot
    sessid = cr.cluster_lookup.sessid(cellid);
    chanid = cr.cluster_lookup.channo(cellid);
    clusonchan = cr.cluster_lookup.clusid(cellid);
    spikes = sessions(sessid).cherries;
    cherryno = find([spikes(:).channr] == chanid & ...
                    [spikes(:).classno] == clusonchan) ;
    assert(numel(cherryno) == 1);
    
    %% sort stimuli according to response-strength in the primed & control condition
    zscore_control = squeeze(zc.zvals(:,:,2));
    zscore_primed = squeeze(zc.zvals(:,:,1));
    [tmemp stimidx] = sort(zscore_control(cellid,:), 'descend');
    [tmemp stimidx_primed] = sort(zscore_primed(cellid,:), 'descend');
    
    %% grab trials of interest & Plot the stimuli
    trials = {};
    conds = [];
    stims = []; 
    snames = {};
    condidx(1,:) = sessions(sessid).condition.condition == 1;
    condidx(2,:) = sessions(sessid).condition.condition == 2;

    for sidx = 1:howmany
        for cidx = 1:2 % condition index
            snames{cidx, sidx} = cr.stim_lookup{stimidx(sidx)};
            tidx = strcmp([snames{cidx,sidx} '.jpg'], sessions(sessid).condition.imagename);
            stims = [stims repmat(stimidx(sidx),1,5)];
            tcidx = tidx & condidx(cidx,:);
            trials = [trials sessions(sessid).cherries(cherryno).trial(tcidx)];
            conds = [conds repmat(cidx, 1, 5)];
            zscorestoplot(sidx, cidx) = zc.zvals(cellid,stimidx(sidx),cidx);
        end
    end
    % write SourceDataFile for ncomms
    if cellid == 3628
        sourceDataFileName = ['source_data_files_ncomms', filesep, 'SourceDataFigure2A.xlsx'];
        
       
        for ti = 1:numel(trials)
            t1.ConditionLabel(ti,1) = conds(ti);
            t1.StimLabel(ti,1) = cr.stim_lookup(stims(ti));
            for ns = 1:numel(trials{ti})
                t1.timestamps(ti,ns) = trials{ti}(ns);
            end
        end
        
        writetable(struct2table(t1), sourceDataFileName, 'Sheet', 'SpikeTimeStamps');
        t2.primedZscores = zscorestoplot(:,1);
        t2.controlZscores = zscorestoplot(:,2);
        t2.StimName = snames(1,:)';
        writetable(struct2table(t2),sourceDataFileName, 'Sheet', 'Zscores');
        
    end

    subplot(1,3,1)
    controlidx = conds == 2;
    plot_raster(trials(controlidx), -1000:2000, stims(controlidx), ...
                repmat([0.2 0.2 0.2; 0 0 0], howmany/2, 1));

    u = 1/(howmany);
    % group trials of same stimulus
    y = [u:u:1-u, u:u:1-u];
    for i = 1:numel(y)
        plot([-1000 2000], [y(i) y(i)], '-k');
    end
    title('control', 'color', 'r');

    %%% do the second row: primed condition
    subplot(1,3,2);
    primedidx = conds == 1;
    plot_raster(trials(primedidx), [-1000:2000]);
    % group trials of same stimulus
    y = [u:u:1-u, u:u:1-u];
    for i = 1:numel(y)
        plot([-1000 2000], [y(i) y(i)], '-k');
    end
    title('primed', 'color', 'b');
    
    %% plot the average zscores next to it
    subplot(1,3,3)
    set(gca, 'color', 'none')
    plot(zscorestoplot(:,1), howmany:-1:1, 'b');
    hold on;
    plot(zscorestoplot(:,2), howmany:-1:1, 'r');
    ylim([0.5 howmany+0.5]);
    xl = xlim;
    for k = 0:howmany
        plot(xl, [k+0.5 k+0.5], '-k');
    end
    set(gca, 'YTick', [1:howmany]);
    set(gca, 'YTickLabel', [howmany:-1:1]);
    xlabel('Z-Scored Firing Rate');
    ylabel('Stimulus #');

    %% plot stimuli next to raster plots
    pos = get(gca, 'Position');
    y = [pos(2):(pos(4)/howmany):(pos(2)+pos(4))];%-pos(4)/howmany;
    hgth = pos(4)/howmany;
    wdth = hgth;
    xp = 0.05;
    for cidx = 1:2
        
        if cidx == 2
           xp = 0.35; 
        end
       
        for i = 1:numel(y)-1
            ax1 = axes('Position', [xp y(i), wdth, hgth]);
            imshow([stimdir snames{cidx, howmany+1-i} '.jpg']);
        end
    end
    set(findall(gcf,'-property','FontSize'),'FontSize',8)
    print(figh, outfn,'-dtiffn','-r600');
    close(figh);
    disp(sprintf('done printing %s',outfn));
end