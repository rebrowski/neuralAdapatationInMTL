function [pos neg] = get_clusters(vals, mode)
%% this is  a subfunction of perm_ttest.m, perm_oneway.m and gets cluster infos
% from a series of values (t/F,1-D array of either from the true or one of the
% false label assignments)
% vals shouls already be thresholded, i.e., only contain vals
% for which p  < clusteralpha, and all other samples should be zero

    if ~exist('mode', 'var')    
        %mode = 'maxsize';   % return only the largest cluster
        mode = 'maxsum'; % return only the cluster with the
                                % highest t-value
        mode = 'all';
    end
    
    %% get positive clusters
    k = vals>0;
    if sum(k)>0
        pos.froms = find(diff(k)>0);
        
        if vals(1)>0
            pos.froms = [1 pos.froms];
        end
        
        pos.tos = find(diff(k)<0);
        
        if vals(end) > 0
            pos.tos = [pos.tos length(vals)];
        end
        pos.sizes = pos.tos - pos.froms;
        for i = 1:length(pos.froms)
            pos.sumts(i) = sum(vals(pos.froms(i):pos.tos(i)));
        end
    else
        pos.froms = [];
        pos.tos = [];
        pos.sizes = [];
        pos.sumts = [];
    end
    
    %% get negative clusters
    k = vals<0;
    if sum(k)>0        
        neg.froms = find(diff(k)>0);
        if vals(1)<0
            neg.froms = [1 neg.froms];
        end
        
        neg.tos = find(diff(k)<0);
     
        if vals(end) < 0
            neg.tos = [neg.tos length(vals)];
        end
        neg.sizes = neg.tos - neg.froms;

        for i = 1:length(neg.froms)
            neg.sumts(i) = sum(vals(neg.froms(i):neg.tos(i)));
        end
    else
        neg.froms = [];
        neg.tos = [];
        neg.sizes = [];
        neg.sumts =[];
    end
    
    posi = [];
    negi = [];
    %% in case only the largest/strongest cluster is desired:
    switch mode
      case 'maxsum';
        posi = find(max(pos.sumts));
        negi = find(min(neg.sumts));
        
      case 'maxsize';
        posi = find(max(pos.sizes));
        negi = find(max(neg.sizes));
    end
    
    if strcmp(mode, 'maxsum') || ...
            strcmp(mode, 'maxsize')
        
        if ~isempty(posi)
            pos.sumts = pos.sumts(posi);
            pos.froms = pos.froms(posi);
            pos.tos = pos.tos(posi);
            pos.sizes = pos.sizes(posi);
        end
        if ~isempty(negi)
            neg.sumts = neg.sumts(negi);
            neg.froms = neg.froms(negi);
            neg.tos = neg.tos(negi);
            neg.sizes = neg.sizes(negi);
        end
    end
    
end