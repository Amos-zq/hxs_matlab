function [ EEG, bcgTemp, bcgTempEpoch ] = bcgRemoval( EEG, etype, method, nc )
%bcgRemoval BCG artifact removal
% Usage:
%   [ EEG, bcgTemp, bcgTempEpoch ] = bcgRemoval( EEG, etype, method, nc )
%
% Inputs:
%   EEG    : EEGLAB data structure
%   etype  : bcg event type
%   mathod : bcg removal method
%   nc     : number of component
% 
% Outputs:
%   EEG          : bcg free EEG structrue
%   bcgTemp      : bcg template of all data
%   bcgTempEpoch : bcg template in epoch
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1: 2013-07-18 17:56, original version 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
bcgEpoch = zeros(EEG.nbchan, 2*PArange+1, length(bcgEvent));
eegEpoch = bcgEpoch;
for i = 1:length(bcgEvent)
    bcgEpoch(:,:,i) = EEGfilted.data(:,bcgEvent(i)-PArange:...
        bcgEvent(i)+PArange);
    eegEpoch(:,:,i) = detrend(EEG.data(:,bcgEvent(i)-PArange:...
        bcgEvent(i)+PArange)', 'constant')';
end
bcgTempEpoch = bcgEpoch;
switch method
    case 'obs-ac'
        pcamat = shiftdim(bcgEpoch,2);
        pcamat = pcamat(:,:)';
        [COEFF, papc] = pca(pcamat);
        papc = papc(:,1:nc);
        for i = 1:length(bcgEvent)
            tempEEG = eegEpoch(:,:,i);
            pad_fit = double(papc)*(double(papc)\double(tempEEG(:)));
            pad_fit = reshape(pad_fit, EEG.nbchan, 2*PArange+1);
            bcgTempEpoch(:,:,i) = reshape(pad_fit, EEG.nbchan, 2*PArange+1);
        end
    case 'tensor'
        % Tensor Decompose BCG
        [C,Z,A] = tensor_BCG(bcgEpoch,nc,50);
        for i = 1:length(bcgEvent)
            bcgTempEpoch(:,:,i) = C*A(:,:,i)*Z;
        end
    case 'sim'
        % SIM Decompose BCG
        [A,S,z] = SIM(bcgEpoch,nc);
        for i = 1:length(bcgEvent)
            bcgTempEpoch(:,:,i) = A*z;
        end
    otherwise
        disp('No such method...');
        return;
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

