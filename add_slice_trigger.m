function EEG = add_slice_trigger( EEG, trType, sliceNum )
%add_slice_trigger Add slice trigger 
%
% Usage:
%   EEG = add_slice_trigger( EEG, trType, sliceNum )
%
% Inputs:
%   EEG:        EEGLAB data structure
%   trType:     TR trigger string
%   sliceNum:   slices in each TR
%
% Outputs:
%   EEG
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   1-Mar-2013, orignal
%   v0.2:   20-Apr-2013, add evaluation of diff tr latency
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

eventTR = struct([]);
eventNum = [];
for i = 1:length(EEG.event)
    if strcmp(EEG.event(i).type, trType)
        eventTR = [eventTR, EEG.event(i)];
        eventNum = [eventn, i];
    end
end
latency = zeros(1, length(eventTR));
for i = 1:length(eventTR)
    latency(i) = eventTR(i).latency;
end
[C,IA,IC] = unique(diff(latency));
if length(C) > 1
    disp('Error: Missing TR')
    for i = 1:length(IA)
        fprintf(1, 'Num = %d, EventNum = %d, TR Length = %d\n', i, eventNum(IA(i)), C(i));
    end
    str = input('Accept? Y/N [Y]: ', 's');
    if strcmp(str, 'N')
        return;
    else
        str = input('Select TR Length Num: ', 's');
        C(1) = C(str2double(str));
    end
end
trLength = C(1);
if rem(trLength, sliceNum)
    disp('Error: Not evenly distributed slices')
    disp('TR:')
    trLength
    disp('sliceNum:')
    sliceNum
    return;
end
sliceInt = trLength/sliceNum;
for i = 1:length(eventTR)
    for j = 1:sliceNum
        sliceEvent = struct('type', 'Slice',...
                            'latency', eventTR(i).latency+sliceInt*(j-1),...
                            'urevent', []);
        EEG.event = [EEG.event, sliceEvent];
    end
end

end

