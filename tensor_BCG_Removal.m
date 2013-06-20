function [ EEG, bcgTemp, C, Z, A, bcgTempEpoch ] = tensor_BCG_Removal( EEG, etype, nc )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% Extract BCG event
bcgEvent = [];
for i = 1:length(EEG.event)
    if strcmp(EEG.event(i).type,etype)
        bcgEvent(end+1) = round(EEG.event(i).latency);
    end
end

% Standard delay between QRS peak and artifact (allen,1998)  
delay = round(0.25*EEG.srate);
bcgEvent = bcgEvent+delay;

% Pulse artifact range
RR=diff(bcgEvent);
mRR=median(RR);
sRR=std(RR);
% PArange=round((0.7*mRR)/2);
% PArange=round((mRR+sRR)/2);
% PArange=round(1.25*(mRR+sRR)/2);
PArange=round(1.5*mRR/2);
midP=PArange+1;

% Delete events out of range
bcgEvent(bcgEvent>(EEG.pnts-PArange)) = [];
bcgEvent(bcgEvent<PArange) = [];

% Filt
EEGfilted = pop_eegfiltnew(EEG, 0.5);
EEGfilted = pop_eegfiltnew(EEGfilted, 0, 40);

% Epoch BCG
bcg = zeros(EEG.nbchan, 2*PArange+1, length(bcgEvent));
for i = 1:length(bcgEvent)
    bcg(:,:,i) = EEGfilted.data(:,bcgEvent(i)-PArange:bcgEvent(i)+PArange);
end

% Tensor Decompose BCG
[C,Z,A] = tensor_BCG(bcg,nc,10);
bcgTempEpoch = bcg;
for i = 1:length(bcgEvent)
    bcgTempEpoch(:,:,i) = C*A(:,:,i)*Z;
end

bcgTemp = zeros(EEG.nbchan, EEG.pnts);
bcgTemp(:,bcgEvent(1)-PArange:bcgEvent(1)+PArange) = bcgTempEpoch(:,:,1);
for i = 2:length(bcgEvent)-1
    PreP = ceil((bcgEvent(i)-bcgEvent(i-1))/2);
    PostP = ceil((bcgEvent(i+1)-bcgEvent(i))/2);
    if PreP > PArange
        PreP = PArange;
    end
    if PostP > PArange
        PostP = PArange;
    end
    bcgTemp(:,bcgEvent(i)-PreP:bcgEvent(i)+PostP) = bcgTempEpoch(:,midP-PreP:midP+PostP,i);
end
bcgTemp(:,bcgEvent(end)-PArange:bcgEvent(end)+PArange) = bcgTempEpoch(:,:,end);

EEG.data = EEG.data - bcgTemp;

end

