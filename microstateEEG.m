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
iterN = 200;
[gfp,gd] = eeg_gfp(EEG.data',1);
peakLoc = peakfinder(zscore(gfp), 1);
[L_gfp, Gamma, alpha_gfp, R_gfp, sigma_mcv_gfp, log] = mscluster((EEG.data(:,peakLoc)), nc, iterN, EEG.chanlocs, 10, 1, 25, 1);

gfp_hrf = mapstd((decimate(conv(double(gfp)', spm_hrf(1/EEG.srate)), EEG.srate*TR, 'FIR')));
gd_hrf = mapstd((decimate(conv(double(gd)', spm_hrf(1/EEG.srate)), EEG.srate*TR, 'FIR')));

[L, alpha] = smoothL(EEG.data, Gamma, iterN, round(EEG.srate*(0.1/2)), 1);

sigma_u = sum((sum(EEG.data.^2, 1)-sum(Gamma(:,L).*EEG.data,1).^2)/(EEG.pnts*(EEG.nbchan-1)), 2);
sigma_mcv = sigma_u*((EEG.nbchan-1)^-1*(EEG.nbchan-1-nc))^-2;
R = sqrt(1 - sigma_u/sum(sum(EEG.data.^2, 1)/(EEG.pnts*(EEG.nbchan-1)), 2));

% conv hrf
alpha_hrf = [];
L_hrf = [];
for i = 1:nc
    alpha_hrf(i,:) = mapstd(decimate(conv(abs(alpha(i,:)), spm_hrf(1/EEG.srate)), EEG.srate*TR, 'FIR'));
    L_hrf(i,:) = mapstd(decimate(conv(double(L==i), spm_hrf(1/EEG.srate)), EEG.srate*TR, 'FIR'));
end


msEEG = struct('gfp', gfp, 'gd', gd, 'gfp_hrf', gfp_hrf, 'gd_hrf', gd_hrf,'peakLoc', peakLoc, ...
               'L', L, 'Gamma', Gamma, 'alpha', alpha, 'R', R, 'sigma_mcv', sigma_mcv, ...
               'alpha_hrf', alpha_hrf, 'L_hrf', L_hrf);

function [L, alpha] = smoothL(V, Gamma, iterN, b, l)
    [Y,L] = max((V'*Gamma).^2,[],2);
    [Ns,Nt] = size(V);
    Nu = size(Gamma,2);
    s_0 = 0;
    eps0 = 10^-6;
    % Initialize progress indicator
    nSteps = 20;
    step = 0;
    fprintf(1, 'relabeling: |');
    strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
    tic
    for n = 1:iterN 
        e = sum((sum(V.^2, 1)-sum((Gamma(:,L).*V).^2, 1))/(Nt*(Ns-1)), 2);
        Nbkt = zeros(Nt, 2*b+1);
        for shift = -b:b
            Nbkt(:,shift+b+1) = circshift(L,shift);
        end
        argk = zeros(Nt,Nu);
        for k = 1:Nu
            argk(:,k) = ((sum(V.^2, 1)-sum((Gamma(:,k)'*V).^2, 1))/(2*e*(Ns-1)))'...
                        - l*sum(Nbkt==k,2);
        end
        [Y,L] = min(argk, [], 2);
        s_u = sum((sum(V.^2, 1)-sum(Gamma(:,L).*V,1).^2)/(Nt*(Ns-1)), 2);
        if abs(s_0-s_u) <= eps*s_u
            fprintf(1,['\n relabel iter stops at ' num2str(n) '\n']);
            break;
        else
            s_0 = s_u;
        end
        [step, strLength] = mywaitbar(n, iterN, step, nSteps, strLength);
    end
    % Deinitialize progress indicator
    fprintf(1, '\n');
    alpha = zeros(Nu, Nt);
    for t = 1:Nt
        kk = L(t);
        for k = 1:Nu
            if k==kk
                alpha(k,t) = V(:,t)'*Gamma(:,k);
            else
                alpha(k,t) = 0;
            end
        end 
    end
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

end

