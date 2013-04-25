function [L, Gamma, alpha, R, sigma_mcv, log] = mscluster(V, Nu, iterN, chanlocs, sm_r, sm_alpha, b, l)
%mscluster EEG microstate clustering analysis
% R. D. Pascual-Marqui, C. M. Michel, and D. Lehmann, ?Segmentation of
% brain electrical activity into microstates: model estimation and
% validation.,? IEEE Trans Biomed Eng, vol. 42, no. 7, pp. 658?665, Jul.
% 1995.
%
% 
% Usage:
%     [L, Gamma, alpha, R, sigma_mcv, log] = mscluster(V, Nu, iterN, chanlocs, sm_r, sm_alpha, b, l);
% 
% Inputs:
%     V:          EEG data, channel*frame
%     Nu:         number of clusters
%     iterN:      number of iteration
%     chanlocs:   chanlocs for spatial smoothness
%     sm_r:       spatial smoothness range
%     sm_alpha:   spatial smoothness ratio
%     b:          time smoothness range
%     l:          time smoothness ratio
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
%     
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% init
[Ns Nt] = size(V);
% 2b)
L = randi(Nu,[Nt 1]);
Gamma = zeros(Ns, Nu);

K = sm_k(chanlocs, sm_r);
% bcg = bcg - repmat(mean(bcg,2), 1, size(bcg,2));

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
%         S = (1+sm_alpha*K)\(V(:,L==k)*V(:,L==k)'/sum(L==k));
        S = (V(:,L==k)*V(:,L==k)'/sum(L==k));
        [X,D] = eig(S);
        [Q,I] = sort(diag(D), 'descend');
        Gamma(:,k) = X(:,I(1));
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

% % relable
% sigma_0 = 0;
% for n = 1:iterN 
%     e = sum((sum(V.^2, 1)-sum((Gamma(:,L).*V).^2, 1))/(Nt*(Ns-1)), 2);
%     Nbkt = zeros(Nt, 2*b+1);
%     for shift = -b:b
%         Nbkt(:,shift+b+1) = circshift(L,b);
%     end
%     argk = zeros(Nt,Nu);
%     for k = 1:Nu
%         argk(:,k) = ((sum(V.^2, 1)-sum((Gamma(:,k)'*V).^2, 1))/(2*e*(Ns-1)))'...
%                     - l*sum(Nbkt==k,2);
%     end
%     [Y,L] = min(argk, [], 2);
%     sigma_u = sum((sum(V.^2, 1)-sum(Gamma(:,L).*V,1).^2)/(Nt*(Ns-1)), 2);
%     if abs(sigma_0-sigma_u) <= eps0*sigma_u
%         fprintf(1,['relabel iter stops at ' num2str(n) '\n']);
%         break;
%     else
%         log(n,2) = abs(sigma_0-sigma_u);
%         sigma_0 = sigma_u;
%     end
% end

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




function K = sm_k(chanlocs, sm_r)
    % spatial smoothness
    nc = size(chanlocs,2);
    chanCor = zeros(6,nc);
    for i = 1:nc
        chanCor(:,i) = [chanlocs(i).Y; chanlocs(i).X; chanlocs(i).Z; chanlocs(i).sph_theta; chanlocs(i).sph_phi; chanlocs(i).sph_radius];
    end
    radius = mean(chanCor(6,:),2);
    G = zeros(nc,nc); Diag = G;
    for i = 1:nc
        for j = 1:nc
    %         G(i,j) = exp(-(1/2)*(sum((chanCor(1:3,i)-chanCor(1:3,j)).^2)/(r^2)));
            G(i,j) = exp(-(1/2)*(distance(chanCor(4,i),chanCor(5,i),chanCor(4,j),chanCor(5,j),radius)^2)/(sm_r^2));
        end
    end
    for i = 1:nc
        Diag(i,i) = sum(G(i,:));
    end;
    K = Diag - G;    
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