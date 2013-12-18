function [ EEG, regCoef, corrER ] = bcgRefLayer( EEG, BCG , chanSel, method, etype)
%BCG removal using reference layer method
% Syntax:  
%     [ EEG, regCoef, corrER ]  = bcgRefLayer( EEG, BCG, method, etype)
%    
%
% Inputs:
%
%
% Outputs:
%
%
% Example:
%     [ EEG, regCoef, corrER ]  = bcgRefLayer( EEG, BCG, method, etype)
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Xiaoshan Huang, xiaoshanhuang@gmail.com
%
% Versions:
%     v0.1: 2013-09-17 16:34, orignal
%     v0.2: 2013-12-06 16:05, add BR, SIM and CCA method
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

chanlocs = BCG.chanlocs;
coorBCG = [extractfield(chanlocs, 'X'); extractfield(chanlocs, 'Y'); extractfield(chanlocs, 'Z');];
chanlocs = EEG.chanlocs;
coorEEG = [extractfield(chanlocs, 'X'); extractfield(chanlocs, 'Y'); extractfield(chanlocs, 'Z');];
dist = zeros(BCG.nbchan, EEG.nbchan);
for i = 1:EEG.nbchan
    dist(:,i) = sqrt(sum((coorBCG - repmat(coorEEG(:,i), [1, size(coorBCG,2)])).^2,1));
end
[Y,sortChan] = sort(dist);
if chanSel > BCG.nbchan
    chanSel = BCG.nbchan;
end
BCG = pop_select(BCG, 'channel', unique(sortChan(1:chanSel,:)));

if nargin < 4
    method = 'AR';
end

if ~strcmp(method, 'AR')
    qrs=[];
    for i=1:length(EEG.event)
        if strcmp(EEG.event(i).type, etype)
            qrs(end+1)=round(EEG.event(i).latency);
        end
    end
    %standard delay between QRS peak and artifact (allen,1998)  
    if strcmp(etype, 'PPG ')
        delay = round(0.4*EEG.srate);
    elseif strcmp(etype, 'qrs')
        delay = round(0.21*EEG.srate);
    end
    mRR = median(diff(qrs));
    qrs(qrs-delay<1) = []; qrs(qrs+mRR-delay>EEG.pnts) = [];
    bcgEpochEEG = epoch(EEG.data, qrs, [-delay mRR-delay]);
    bcgEpochBCG = epoch(BCG.data, qrs, [-delay mRR-delay]);
end
corrER = [];
regCoef = zeros(BCG.nbchan, EEG.nbchan);
switch method
    case 'AR'
        for chan = 1:EEG.nbchan
            regCoef(:, chan) = BCG.data'\EEG.data(chan,:)';
            refSignal = regCoef(:, chan)' * BCG.data;
            corrER(chan) = corr(EEG.data(chan,:)', refSignal');
            EEG.data(chan,:) = EEG.data(chan,:) - refSignal;
        end
    case 'BR'
        for chan = 1:EEG.nbchan
            regCoef(:, chan) = bcgMeanBCG'\bcgMeanEEG(chan,:)';
            refSignal = regCoef(:, chan)' * BCG.data;
            corrER(chan) = corr(EEG.data(chan,:)', refSignal');
            EEG.data(chan,:) = EEG.data(chan,:) - refSignal;
        end
    case 'SIMR'
        [simA, simS, simz] = SIM(bcgEpochBCG);
        comps = 1:20;
        for chan = 1:EEG.nbchan
            regCoef(:, chan) = (simA(:,comps)*simS(comps,:)*bcgEpochBCG(:,:))' \ bcgEpochEEG(chan,:)';
            refSignal = regCoef(:, chan)' * simA(:,comps)*simS(comps,:)*BCG.data;
            corrER(chan) = corr(EEG.data(chan,:)', refSignal');
            EEG.data(chan,:) = EEG.data(chan,:) - refSignal;
        end
    case 'CCAR'
        bcgMeanEEG = mean(bcgEpochEEG, 3);
        bcgMeanBCG = mean(bcgEpochBCG, 3);
        numComp = min(EEG.nbchan, BCG.nbchan);
        [coeff, pcs] = pca(bcgMeanBCG', 'NumComponents', numComp);
        [ccaSBCG,ccaSEEG,ccaR,U,V] = canoncorr(pcs, bcgMeanEEG');
        ccaSBCG = coeff * ccaSBCG;
        ccaABCG = pinv(ccaSBCG);
        ccaAEEG = pinv(ccaSEEG);
        comps = find(ccaR > 0.9);
        for chan = 1:EEG.nbchan
            regCoef(:, chan) = ccaSBCG(:,comps)*ccaAEEG(comps,chan);
            refSignal = regCoef(:, chan)' * BCG.data;
            corrER(chan) = corr(EEG.data(chan,:)', refSignal');
            EEG.data(chan,:) = EEG.data(chan,:) - refSignal;
        end
    case 'CCAMW'
        winSize = 30 * EEG.srate;
        stepLength = 1 * EEG.srate;
        for pnt = 1:stepLength:EEG.pnts
            if pnt <= winSize
                bcgEpochRange = find(qrs<winSize);
            else
                bcgEpochRange = find(qrs>pnt-winSize & qrs<pnt);
            end
            bcgMeanEEG = mean(bcgEpochEEG(:,:,bcgEpochRange), 3);
            bcgMeanBCG = mean(bcgEpochBCG(:,:,bcgEpochRange), 3);
%             bcgMeanEEG = rSIM(bcgEpochEEG(:,:,bcgEpochRange));
%             bcgMeanBCG = rSIM(bcgEpochBCG(:,:,bcgEpochRange));
            numComp = min(EEG.nbchan, BCG.nbchan);
            [coeff, pcs] = pca(bcgMeanBCG', 'NumComponents', numComp);
            [ccaSBCG,ccaSEEG,ccaR,U,V] = canoncorr(pcs, bcgMeanEEG');
            ccaSBCG = coeff * ccaSBCG;
            ccaAEEG = pinv(ccaSEEG);
            comps = find(ccaR > 0.9);
            regCoef = ccaSBCG(:,comps)*ccaAEEG(comps,:);
            if pnt+stepLength-1 > EEG.pnts
                pntRange = pnt:EEG.pnts;
            else
                pntRange = pnt:pnt+stepLength-1;
            end
            EEG.data(:,pntRange) = EEG.data(:,pntRange) - regCoef' * BCG.data(:,pntRange);
        end
    case 'RBCCA'
        winSize = 30 * EEG.srate;
        stepLength = 1 * EEG.srate;
        for pnt = 1:stepLength:EEG.pnts
            if pnt <= winSize
                bcgTempRange = 1:winSize;
            else
                bcgTempRange = pnt-winSize:pnt ;
            end
            [coeff, bcgPCs] = pca(BCG.data(:,bcgTempRange)');
            pcRange = 1:min(BCG.nbchan, EEG.nbchan);
            coeff = coeff(:,pcRange);
            bcgPCs = bcgPCs(:,pcRange)';
%             [ccaSBCG,ccaSEEG,ccaR] = canoncorr(bcgPCs', EEG.data(:,bcgTempRange)');
            [ccaSBCG,ccaSEEG,ccaR] = rbcca(bcgPCs, EEG.data(:,bcgTempRange));
            ccaSBCG = coeff(:,pcRange) * ccaSBCG;
            ccaAEEG = pinv(ccaSEEG);
            comps = find(ccaR > 0.9);
            regCoef = ccaSBCG(:,comps)*ccaAEEG(comps,:);
            if pnt+stepLength-1 > EEG.pnts
                pntRange = pnt:EEG.pnts;
            else
                pntRange = pnt:pnt+stepLength-1;
            end
            EEG.data(:,pntRange) = EEG.data(:,pntRange) - regCoef' * BCG.data(:,pntRange);
        end
end
EEG.setname = [EEG.setname '_BRL_' method];

