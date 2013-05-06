function msEEG = microstateEEG( EEG, nc )
% microstateEEG calculate EEG microstates using mscluster cluster
% Usage:
%   msEEG = microstateEEG( EEG, nc );
% 
% Inputs:
%   EEG :EEGLAB data structure
%   nc  :Number of Clusters
% 
% Outputs:
%   msEEG.
%   gfp     :gfp
%   gd      :gd
%   peakLoc :gfp peak location
%   L       :cluster labels
%   Gamma   :cluster map
%   alpha   :time course of each cluster
%   
% 
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
% 
% Versions:
%   v0.1:   27-Mar-2013, orignal
%   v0.2:   29-Mar-2013, add hrf conv
%   v0.3:   06-May-2013, use mscluster
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TR = 2.04;
[gfp,gd] = eeg_gfp(EEG.data',1);
peakLoc = peakfinder(zscore(gfp), 1);
[L_gfp, Gamma, alpha_gfp, R_gfp, sigma_mcv_gfp, log] = mscluster(mapstd(EEG.data(:,peakLoc)), nc, 200, EEG.chanlocs, 10, 1, 25, 1);

gfp_hrf = mapstd((decimate(conv(double(gfp)', spm_hrf(1/EEG.srate)), EEG.srate*TR, 'FIR')));
gd_hrf = mapstd((decimate(conv(double(gd)', spm_hrf(1/EEG.srate)), EEG.srate*TR, 'FIR')));

% labeling all data points, calculate alpha for all data points
[Y,L] = max((EEG.data'*Gamma).^2,[],2);
alpha = zeros(nc, EEG.pnts);
for t = 1:EEG.pnts
    kk = L(t);
    for k = 1:nc
        if k==kk
            alpha(k,t) = EEG.data(:,t)'*Gamma(:,k);
        else
            alpha(k,t) = 0;
        end
    end 
end
sigma_u = sum((sum(EEG.data.^2, 1)-sum(Gamma(:,L).*EEG.data,1).^2)/(EEG.pnts*(EEG.nbchan-1)), 2);
sigma_mcv = sigma_u*((EEG.nbchan-1)^-1*(EEG.nbchan-1-nc))^-2;
R = sqrt(1 - sigma_u/sum(sum(EEG.data.^2, 1)/(EEG.pnts*(EEG.nbchan-1)), 2));

% conv hrf
alpha_hrf = [];
L_hrf = [];
for i = 1:nc
    alpha_hrf(i,:) = decimate(conv(abs(alpha(i,:)), spm_hrf(1/EEG.srate)), EEG.srate*TR, 'FIR');
    L_hrf(i,:) = decimate(conv(double(L==i), spm_hrf(1/EEG.srate)), EEG.srate*TR, 'FIR');
end


msEEG = struct('gfp', gfp, 'gd', gd, 'gfp_hrf', gfp_hrf, 'gd_hrf', gd_hrf,'peakLoc', peakLoc, ...
               'L', L, 'Gamma', Gamma, 'alpha', alpha, 'R', R, 'sigma_mcv', sigma_mcv, ...
               'alpha_hrf', alpha_hrf, 'L_hrf', L_hrf);

end

