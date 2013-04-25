function [ pc, times] = eeg_pca( EEG, events, tlim, npc )
%eeg_pca perform pca on EEG event epochs
% Usage:
%   [ pc times] = eeg_pca( EEG, events, tlim, npc )
%
% Inputs:
%   EEG:        EEGLAB data structure
%   events:     event type string
%   tlim:       time limits for epoch, in s
%   npc:        number of principle components
% Outputs:
%   pc:         pc
%   times:      times of epoch
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   04-Apr-2013, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

bcgEEG = pop_epoch( EEG, {events}, tlim);
pc = zeros(bcgEEG.nbchan, bcgEEG.pnts, npc);
for i = 1:bcgEEG.nbchan
    pcamat = squeeze(bcgEEG.data(i,:,:))';
    pcs = npc-1;
    pcamat=detrend(pcamat','constant')';
    meaneffect=mean(pcamat);
    dpcamat=detrend(pcamat,'constant');
    [apc,ascore,asvar]=pca_calc(dpcamat');
    papc=[meaneffect' ascore(:,1:pcs)];
    pc(i,:,:) = papc;
end
times = bcgEEG.times;
figure
for i = 1:npc
    subplot(npc,1,i), plot(times, pc(:,:,i)), title(['PC' num2str(i)]);
end
end

