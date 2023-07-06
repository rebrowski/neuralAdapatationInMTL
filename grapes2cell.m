function [outspikes] = grapes2cell(inspikes)
%% function [spikes] = convert_spikearray2struct(spikes1)
% takes an array of spike-times stored in grapes.mat (eg. chan42.class2.stim{4})
% removes the funny 10000 stores spike times in a cell

ntrials = size(inspikes, 1);
for t = 1:ntrials
    sp = inspikes(t,inspikes(t,:)~= 10000);
    
    % this should eliminate some of the artefacts
    if numel(sp) > 200
        sp = ones(0,1);
    end
    
    outspikes{t} = sp;
end

if ~exist('outspikes','var')
    outspikes = {};
end
