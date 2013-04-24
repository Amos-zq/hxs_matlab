function msEEG = microstateEEG( EEG, nc )
%microstateEEG calculate EEG microstates using kmeans cluster
% Usage:
%   msEEG = microstateEEG( EEG, nc );
%
% Inputs:
%   EEG:        EEGLAB data structure
%   nc:         Number of Clusters
%
% Outputs:
%   msEEG.
%   gfp:        gfp
%   gd:         gd
%   peakLoc:    gfp peak location
%   microstate: microstate
%   IDX:        cluster index
%   C:          cluster map
%   msCorr:     corr with data
%   msCluster:  cluster index in original dimension
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   27-Mar-2013, orignal
%   v0.2:   29-Mar-2013, add hrf conv
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[gfp,gd] = eeg_gfp(EEG.data',1);
peakLoc = peakfinder(zscore(gd), 1);
microstate = mapstd(EEG.data(:,peakLoc)')';
fprintf(1, 'Clustering...\n');
[IDX, C] = kmeans(microstate', nc);
R = corrcoef(C');
[I, J] = find((tril(R)<-0.9)==1);
while ~isempty(I)
    for i = 1:length(I)
        microstate(:, IDX==I(i)) = -microstate(:, IDX==I(i));
    end
    fprintf(1, 'Clustering...\n');
    [IDX, C] = kmeans(microstate', nc);
    R = corrcoef(C');
    [I, J] = find((tril(R)<-0.99)==1);
end

gfp_hrf = conv(mapstd(decimate(double(gfp'), 250*2.04)), spm_hrf(2.04));
gd_hrf = conv(mapstd(decimate(double(gd'), 250*2.04)), spm_hrf(2.04));


msCluster = zeros(nc, EEG.pnts);
for i = 1:nc
    msCluster(i,peakLoc(IDX==i)) = 1;
end

msCluster_hrf = [];
for i = 1:nc
    msCluster_hrf(i,:) = conv(mapstd(decimate(double(msCluster(i,:)), 250*2.04)), spm_hrf(2.04)); 
end

msCorr = zeros(nc, EEG.pnts);
msDiss = msCorr;

% Initialize progress indicator
nSteps = 20;
step = 0;
fprintf(1, 'corr: |');
strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
tic
for i = 1:EEG.pnts
    for j = 1:nc
        msCorr(j,i) = corr(EEG.data(:,i), C(j,:)'); 
        msDiss(j,i) = sqrt(2*(1-msCorr(j,i)));
    end
    [step, strLength] = mywaitbar(i, EEG.pnts, step, nSteps, strLength);
end
% Deinitialize progress indicator
fprintf(1, '\n')

msCorr_hrf = [];
msDiss_hrf = [];
for i = 1:nc
    msCorr_hrf(i,:) = conv(mapstd(decimate(abs(msCorr(i,:)), 250*2.04)), spm_hrf(2.04));
    msDiss_hrf(i,:) = conv(mapstd(decimate(abs(msDiss(i,:)), 250*2.04)), spm_hrf(2.04));
end

msEEG = struct('gfp', gfp, 'gd', gd, 'gfp_hrf', gfp_hrf, 'gd_hrf', gd_hrf, ...
               'peakLoc', peakLoc, 'microstate', microstate, ...
               'IDX', IDX, 'C', C, 'msCluster', msCluster, 'msCluster_hrf', msCluster_hrf, ...
               'msCorr', msCorr, 'msCorr_hrf', msCorr_hrf, 'msDiss', msDiss, 'msDiss_hrf', msDiss_hrf);

end




function [step, strLength] = mywaitbar(compl, total, step, nSteps, strLength)

progStrArray = '/-\|';
tmp = floor(compl / total * nSteps);
if tmp > step
    fprintf(1, [repmat('\b', 1, strLength) '%s'], repmat('=', 1, tmp - step))
    step = tmp;
    ete = ceil(toc / step * (nSteps - step));
    strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '%s %3d%%, ETE %02d:%02d'], progStrArray(mod(step - 1, 4) + 1), floor(step * 100 / nSteps), floor(ete / 60), mod(ete, 60));
end

end

