function [EEG, fitted_art, papc_ac] = obs_ac( EEG,etype,npc )
%pulse artifacts removal using Optimal Basis Set from ALL Channel
%
% Usage:
%   [EEG, fitted_art, papc_ac] = obs_ac( EEG,etype,npc )
%
% Inputs:
%   EEG:    EEGLAB data structure
%   etype:  bcg event type, e.g. 'bcg'
%   npc:    number of OBS to remove
%
% Outputs:
%   EEG:    clear EEG
%   fitted_art: fitted artifact
%   papc_ac:    artifact PCs from ALL Channel
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   19-Jun-2013, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%init
%-----
fs=EEG.srate; 
[channels, samples]=size(EEG.data);
pcs=npc-1;

%standard delay between QRS peak and artifact (allen,1998)  
delay=round(0.21*fs);

Gwindow=20; 
GHW=floor(Gwindow/2);
rcount=0;
firstplot=1;

% memory allocation
%------------------
avgdata=zeros(Gwindow,round(fs*2));
drift=zeros(1,round(fs*2)*Gwindow);
fitted_art=zeros(channels,samples);
mPA=zeros(1,round(2*fs));
peakplot=zeros(1,samples);
mEye=eye(channels);

%Extract QRS events
%------------------

QRSevents=[];
for E=1:length(EEG.event)
    if strcmp(EEG.event(E).type,etype)
        QRSevents(end+1)=round(EEG.event(E).latency);
    end
end

peakplot(QRSevents)=1;
sh=zeros(1,delay);
np=length(peakplot);
peakplot=[sh peakplot(1:np-delay)];
peakI=find(peakplot>0);
peakcount=length(peakI);
pa=peakcount;

RR=diff(peakI);
mRR=median(RR);
sRR=std(RR);
PArange=round((0.7*mRR)/2);
% PArange=round((mRR+sRR)/2);
% PArange=round(1.25*(mRR+sRR)/2);
% PArange=round(1.5*mRR/2);
midP=PArange+1;
while (peakI(pa)+PArange > samples)
    pa=pa-1;
end
steps=channels*pa;
peakcount=pa;

%make filter
%-----------
a=[0 0 1 1];
f=[0 0.4/(fs/2) 0.9/(fs/2) 1];
ord=round(3*fs/0.5);
fwts=firls(ord,f,a);

%Channel-wise obs
eegfilted = zeros(channels, samples);
for ch=1:channels
    eegfilted(ch,:)=filtfilt(fwts,1,double(EEG.data(ch,:)));
end
pcamat=zeros(channels,2*PArange+1,pa-1);
for p=2:pa
    pcamat(:,:,p-1)=eegfilted(:,peakI(p)-PArange:peakI(p)+PArange);
end
pcamat = shiftdim(pcamat,2);
pcamat = pcamat(:,:);
pcamat=detrend(pcamat','constant')';
meaneffect=mean(pcamat);
dpcamat=detrend(pcamat,'constant');
[apc,ascore,asvar]=pca_calc(dpcamat');
papc_ac = [meaneffect' ascore(:,1:pcs)];
papc_ac = reshape(papc_ac, channels, 2*PArange+1, npc);

%Artifact Subtraction
%---------------------

% for the first window/2 points use arthemitic mean for averageing.
% findg mean QRS peak-to-peak (R-to-R) interval
for ch=1:channels
    
    papc = squeeze(papc_ac(ch, :, :));
    
    try
        pad_fit=double(papc)*(double(papc)\...
            double(detrend(EEG.data(ch,peakI(1)-PArange:...
            peakI(1)+PArange)','constant')));

        fitted_art(ch,peakI(1)-PArange:peakI(1)+PArange)=...
            pad_fit';
    catch
    end

    %-----------------------------------------------------
    for p=2:GHW+1
       pad_fit=double(papc)*(double(papc)\...
            double(detrend(EEG.data(ch,peakI(p)-PArange:...
            peakI(p)+PArange)','constant')));

       fitted_art(ch,peakI(p)-PArange:peakI(p)+PArange)=...
            pad_fit';
    end

    %---------------- Processing of central data ---------------------
    %cycle through peak artifacts identified by peakplot
    rcount=GHW;
    for p=GHW+2:peakcount-GHW-1 
        PreP=ceil((peakI(p)-peakI(p-1))/2);
        PostP=ceil((peakI(p+1)-peakI(p))/2);
        if PreP > PArange
            PreP=PArange;
        end
        if PostP > PArange
            PostP=PArange;
        end

       pad_fit=double(papc(midP-PArange:midP+PArange,:))*...
           (double(papc(midP-PArange:midP+PArange,:))\...
            double(detrend(EEG.data(ch,peakI(p)-PArange:...
            peakI(p)+PArange)','constant')));

       fitted_art(ch,peakI(p)-PreP:peakI(p)+PostP)=...
            pad_fit(midP-PreP:midP+PostP)';
    end

    %---------------- Processing of last GHW peaks------------------
    sectionPoints=samples-(peakI(peakcount-GHW)-PreP)+1;

    for p=peakcount-GHW:peakcount
        try

           pad_fit=double(papc(midP-PArange:midP+PArange,:))*...
               (double(papc(midP-PArange:midP+PArange,:))\...
                double(detrend(EEG.data(ch,peakI(p)-PArange:...
                peakI(p)+PArange)','constant')));

           fitted_art(ch,peakI(p)-PreP:peakI(p)+PostP)=...
                pad_fit(midP-PreP:midP+PostP)';


        catch
        end        
    end

end

EEG.data=EEG.data-fitted_art;

end

