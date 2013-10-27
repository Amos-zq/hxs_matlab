function [L, Gamma, alpha, R, sigma_mcv, log] = mscluster(V, Nu, iterN, ts, ss, bcg)
%mscluster EEG microstate clustering analysis
% R. D. Pascual-Marqui, C. M. Michel, and D. Lehmann, ?Segmentation of
% brain electrical activity into microstates: model estimation and
% validation.,? IEEE Trans Biomed Eng, vol. 42, no. 7, pp. 658?665, Jul.
% 1995.
%
% 
% Usage:
%     [L, Gamma, alpha, R, sigma_mcv, log] = mscluster(V, Nu, iterN, ts, ss, bcg);
% 
% Inputs:
%     V:          EEG data, channel*frame
%     Nu:         number of clusters
%     iterN:      number of iteration
%     ts:         time smoothness range and ratio [ts_b ts_l]
%     ss:         spatial smoothness kernal ss_alpha*ss_K
%     bcg:        bcg noise
% 
% Outputs:
%     L:          microstate label
%     Gamma:      microstate spatial map
%     alpha:      microstate time course
%     R:          explained variance
%     sigma_mcv:  cross-validation residual variance
%     log:        interation log
%    
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
% 
% Versions:
%     v0.1:   01-Apr-2013, original
%     v0.2:   24-Apr-2013, revesion
%     v0.3:   21-May-2013, optimize relabeling
%     v0.4:   22-May-2013, add bcg
%     v0.5:   01-Jul-2013, optimize input params
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 4
    ts = []; ss = []; bcg = [];
end
if nargin < 5
    ss = []; bcg = [];
end
if nargin < 6
    bcg = [];
end

if isempty(ts)
    time_smooth = 0;
else
    time_smooth = 1;
    ts_b = ts(1);
    ts_l = ts(2);
end

if isempty(ss)
    spatial_smooth = 0;
else
    spatial_smooth = 1;
    ss_K = ss;
end
if isempty(bcg)
    bcg_reduction = 0;
else
    bcg_reduction = 1;
    bcg = bcg - repmat(mean(bcg,2), 1, size(bcg,2));
    ECovRes = bcg*bcg'/size(bcg,2);
end

% init
[Ns,Nt] = size(V);
% 2b)
L = randi(Nu,[Nt 1]);
Gamma = zeros(Ns, Nu);
eps0 = 10^-6;
sigma_0 = 0; 
% Initialize progress indicator
nSteps = 20;
step = 0;
fprintf(1, 'clustering: |');
strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
tic
for n = 1:iterN 
    % Calculate Gamma
    for k = 1:Nu
        S = (V(:,L==k)*V(:,L==k)'/sum(L==k));
        if spatial_smooth
            S = (1+ss_K)\(V(:,L==k)*V(:,L==k)'/sum(L==k));
        end
        if bcg_reduction
            Wh = ECovRes^(-0.5);
            S = (Wh*V(:,L==k))*(Wh*V(:,L==k))'/sum(L==k);
        end
        
        [X,D] = eig(S);
        [Q,I] = sort(diag(D), 'descend');
        Gamma(:,k) = X(:,I(1));
        if bcg_reduction
            Gamma(:,k) = pinv(A(:,k));
        end
    end
    sigma_u = sum((sum(V.^2, 1)-sum(Gamma(:,L).*V,1).^2)/(Nt*(Ns-1)), 2);
    if abs(sigma_0-sigma_u) <= eps*sigma_u
        fprintf(1,['\n cluster iter stops at ' num2str(n) '\n'])
        break;
    else
        log(n,1) = abs(sigma_0-sigma_u);
        sigma_0 = sigma_u;
        [Y,L] = max((V'*Gamma).^2,[],2);
    end
    [step, strLength] = mywaitbar(n, iterN, step, nSteps, strLength);
end
% Deinitialize progress indicator
fprintf(1, '\n');

if time_smooth
    % relable
    sigma_0 = 0;
    % Initialize progress indicator
    nSteps = 20;
    step = 0;
    fprintf(1, 'relabeling: |');
    strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
    tic
    for n = 1:iterN 
        e = sum((sum(V.^2, 1)-sum((Gamma(:,L).*V).^2, 1))/(Nt*(Ns-1)), 2);
        Nbkt = zeros(Nt, 2*ts_b+1);
        for shift = -ts_b:ts_b
            Nbkt(:,shift+ts_b+1) = circshift(L,shift);
        end
        argk = zeros(Nt,Nu);
        for k = 1:Nu
            argk(:,k) = ((sum(V.^2, 1)-sum((Gamma(:,k)'*V).^2, 1))/(2*e*(Ns-1)))'...
                        - ts_l*sum(Nbkt==k,2);
        end
        [Y,L] = min(argk, [], 2);
        sigma_u = sum((sum(V.^2, 1)-sum(Gamma(:,L).*V,1).^2)/(Nt*(Ns-1)), 2);
        if abs(sigma_0-sigma_u) <= eps0*sigma_u
            fprintf(1,['\n relabel iter stops at ' num2str(n) '\n']);
            break;
        else
            log(n,2) = abs(sigma_0-sigma_u);
            sigma_0 = sigma_u;
        end
        [step, strLength] = mywaitbar(n, iterN, step, nSteps, strLength);
    end
    % Deinitialize progress indicator
    fprintf(1, '\n');
end

sigma_mcv = sigma_u*((Ns-1)^-1*(Ns-1-Nu))^-2;

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

R = 1 - sigma_u/sum(sum(V.^2, 1)/(Nt*(Ns-1)), 2);
R = sqrt(R);


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