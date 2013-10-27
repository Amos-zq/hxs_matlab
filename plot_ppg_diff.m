function [ ppgLatency, ppgIndex, rri ] = plot_ppg_diff (EEG, ppgType, rriVar)
%plot_ppg_diff plot time interval between ppg trigger
%
% Usage:
%   [ ppgIndex, ppgDiff, rri ] = plot_ppg_diff (EEG, ppgType, rriVar)
%
% Inputs:
%   EEG:        EEGLAB data structure
%   ppgType:    PPG trigger string
%   rriVar:     RR interval variation range [low high]
% Outputs:
%   ppgLatency: ppg event latency
%   ppgIndex:   ppg event index
%   ppgDiff:    time interval between ppg triggers
%   rri:        estimated rr interval
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   24-Mar-2013, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ppgLatency = []; ppgIndex = [];
for i = 1:length(EEG.event)
    if strcmp(EEG.event(i).type, ppgType)
        ppgLatency = [ppgLatency EEG.event(i).latency];
        ppgIndex = [ppgIndex i];
    end
end
ppgDiff = diff(ppgLatency)';
z = linkage(ppgDiff);
rri = mean(ppgDiff(z(z(:,3)==0, 1)));

ppgDiffPlot = [ppgDiff ...
               rri*ones(length(ppgDiff),1) ...
               (1-rriVar(1))*rri*ones(length(ppgDiff),1) ...
               (1+rriVar(2))*rri*ones(length(ppgDiff),1)]/EEG.srate;
figure, plot(ppgDiffPlot);

end