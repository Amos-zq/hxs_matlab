function glm = myglm( X, func, mask, alpha, ifr, id )
%myglm glm analysis of nii img, given a design matrix
%
% Usage:
%   glm = myglm( X, func, mask, ifr, id )
%
% Inputs:
%   X:      design matrix
%   func:   func  
%   alpha:  alpha threshold
%   ifr:    frame index to analysis, [start end]
%   id:     design index to analysis
%
% Outputs:
%   glm:
%   b:  beta
%   t:  tscore
%   p:  p value, p_fdr, p_masked
%   r:  res
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   12-Mar-2013, orignal
%   v0.2:   29-Mar-2013, replace nii with func, add mask
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isempty(ifr)
    ifr(1) = 1;
    ifr(2) = size(func, 4);
end
nfr = ifr(2) - ifr(1) +1;
% nii = load_nii(fpnii);
% X = load(fpX);
[nf nd] = size(X);

% option to use some of design
if (~isempty(id))
    X = X(:, id);
end

fprintf('Using %d convariates\n', size(X,2));

% add second order detrending covs

c1 = (1:nf)';		% detrend with quadratic
c2 = c1.*c1;
X = [X c1 c2];

X = X(ifr(1):ifr(2),:);

[nx, ny, nslc, nf] = size(func);
if isempty(mask)
    mask = ones(nx, ny, nslc);
end


% Initialize progress indicator
maskIdx = find(mask==1);
nTotal = length(maskIdx);
nSteps = 20;
step = 0;
fprintf(1, 'myglm(): |');
strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
tic
b = zeros(nx, ny, nslc, size(X,2)+1);
t = b;
p = b;
r = zeros(nx, ny, nslc);

% begin glm loop
for islc = 1:nslc
    for ix = 1:nx
        for iy = 1:ny
            if mask(ix,iy,islc)
                Y = double(squeeze(func(ix,iy,islc,ifr(1):ifr(2))));
                [B, dev, Stat] = glmfit(X, Y, 'normal');
                b(ix, iy, islc, :) = B;
                Stat.t(isnan(Stat.t)) = 0;
                t(ix, iy, islc, :) = Stat.t;
                p(ix, iy, islc, :) = Stat.p;
                r(ix, iy, islc) = sqrt(dev);
                [step, strLength] = mywaitbar(find(maskIdx==sub2ind([nx,ny,nslc],ix,iy,islc)), nTotal, step, nSteps, strLength);
            end
        end
    end
end

% Deinitialize progress indicator
fprintf(1, '\n')

if isempty(alpha)
    p_fdr = []; p_masked = [];
else
    [p_fdr, p_masked] = fdr(p, alpha);
    t = t.*p_masked;
    p = p.*p_masked;
end
glm = struct('b', b, 't', t, 'p', p, 'p_fdr', p_fdr, 'p_masked', p_masked, 'r', r);

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
