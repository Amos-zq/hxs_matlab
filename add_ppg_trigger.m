function [ EEGOUT ] = add_ppg_trigger( EEG, ppgType, rriVar )
%add_ppg_trigger add missing ppg trigger
% Usage:
%   EEG = add_ppg_trigger( EEG, 'PPG ' )
%
% Inputs:
%   EEG:        EEGLAB data structure
%   ppgType:    PPG trigger string
%   rriVar:     RR interval variation range [low high]
% Outputs:
%   EEG
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   1-Mar-2013, orignal
%   V0.2:   22-Mar-2013, add plot orignal ppg latency
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

EEG = pop_editeventvals( EEG, 'sort', { 'latency', 0 } );
[ ppgLatency, ppgIndex, rri ] = plot_ppg_diff (EEG, ppgType, rriVar);
pause;

for i = 2:length(ppgLatency)-1
    il = ppgLatency(i) - ppgLatency(i-1);
    ir = ppgLatency(i+1) - ppgLatency(i);
    rl = (1-rriVar(1))*rri;
    rh = (1+rriVar(2))*rri;
    if (il<rl || il>rh) || ((ir<rl || ir>rh))
        EEG.event(ppgIndex(i)).type = 'bPPG';
    end
end

[ ppgLatency, ppgIndex, rri ] = plot_ppg_diff (EEG, ppgType, rriVar);
pause;

for i = 1:length(ppgLatency)-1
    idiff = ppgLatency(i+1) - ppgLatency(i);
    nrr = round(idiff/rri);
    rri_temp = round(idiff/nrr);
    if nrr >1 
        for j = 1:(nrr-1)
            event_add = EEG.event(ppgIndex(i));
            event_add.latency = event_add.latency+j*rri_temp;
            EEG.event = [EEG.event event_add];
        end
    end
end

EEG = pop_editeventvals( EEG, 'sort', { 'latency', 0 } );
[ ppgLatency, ppgIndex, rri ] = plot_ppg_diff (EEG, ppgType, rriVar);
pause;

EEGOUT = EEG;

end

