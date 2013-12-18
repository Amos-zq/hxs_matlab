function [erp, snrMean, peakAmp] = erpsnr( EEG, eventTypes, epochLimits, peakRange )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
%
% Syntax:  
%     
%
% Inputs:
%     
%
% Outputs:
%     
%
% Example:
%     
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Xiaoshan Huang, xiaoshanhuang@gmail.com
%
% Versions:
%    v0.1: , orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbPeak = size(peakRange,1);
nbEvent = length(eventTypes);
idxPeak = round((peakRange - epochLimits(1))/1000*EEG.srate);
epochEEG = cell(nbEvent,1);

for event = 1:nbEvent
    epoch = pop_epoch(EEG, {eventTypes{event}}, epochLimits/1000);
    epochEEG{event} = pop_rmbase(epoch, [epochLimits(1) 0]);
end

nbchan = EEG.nbchan;
snrMean = zeros(nbchan, nbPeak, nbEvent);
pnts = epochEEG{1}.pnts;
erp = zeros(nbchan, pnts, nbEvent);
peakAmp = cell(nbEvent,1);

for event = 1:nbEvent
    epoch = epochEEG{event};
    erp(:,:,event) = mean(epoch.data,3);
    erpNrom = erp(:,:,event) ./ repmat(std(erp(:,epoch.times<0,event),0,2), [1,pnts]);
    peakAmpTrial = zeros(nbchan, nbPeak, epoch.trials);
    for peak = 1:nbPeak
        index = idxPeak(peak,1):idxPeak(peak,2);
        snrMean(:,peak,event) = 10*log10(mean(erpNrom(:,index).^2,2));
        trialNorm = epoch.data ./ repmat(std(epoch.data(:,epoch.times<0,:),0,2), [1,pnts,1]);
        peakAmpTrial(:,peak,:) = mean(trialNorm(:,index,:),2);
    end
    peakAmp{event} = peakAmpTrial;
end


end

