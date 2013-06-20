function [C,Z,A] = tensor_BCGremove(BCG,nc,iter)
% input : 
% BCG nc iter
% BCG = chan*pecelen*trial each should be zero mean
% nc: number of components
% iter: iteration times
% output :
% C spatial pattern: chan*components 
% Z BCG signal: components*pacelen
% A amplitude matrix: components*components*trial, the elements on diagonal
% is the amplitude of components
 
% my model:Xi = C*Z*ai+e
% |C|=1
% |Z|=1
% e~N(0,phi)


%%%% initialize %%%%%%%%%%%%%%%%
if (nargin<1||nargin>3)
    error('input error');
  
end
if nargin == 1
    nc =4;
    iter =200;
end
if nargin == 2
    iter = 200;
end

[chan,pacelen,trial] = size(BCG);
% phi
phi = zeros(chan);
BCGmean = mean(BCG,3);
for ipeak = 1:trial
    temp = BCG(:,:,ipeak)-BCGmean;
    phi = phi+temp*temp';
end
phi = phi./double(trial*pacelen);
% ai
A = zeros(nc,nc,trial);

for ipeak = 1:trial
    A(:,:,ipeak) = eye(nc);
end
% C,Z
[U,S,V] = svd(BCGmean,'econ');
C = U(:,1:nc);
Z = S(1:nc,1:nc)*V(:,1:nc)';

for it = 1:iter
%%%%% iteration %%%%%%%%%%%%%%%
phi_wh = phi^-0.5;
% X_wh = phi_wh*BCG;
% % P = CP_ALS(X_wh,nc);
% 
% C_wh = phi_wh*C;
%%%% update C %%%%%%%
tempU = (A(:,:,1)*Z)*(A(:,:,1)*Z)';
tempW = BCG(:,:,1)*(A(:,:,1)*Z)';
for ipeak = 2:trial
    tempU = tempU + (A(:,:,ipeak)*Z)*(A(:,:,ipeak)*Z)';
    tempW = tempW + BCG(:,:,ipeak)*(A(:,:,ipeak)*Z)';
end
C = tempW*inv(tempU);

%%%% update A %%%%%
C_wh = phi_wh*C;
% ttm
% X_wh = phi_wh*BCG;
C_inv = pinv(C_wh);
for ipeak = 1:trial
    X_wh = phi_wh*BCG(:,:,ipeak);
    tempX = C_inv * X_wh;
    for ic = 1: nc
       
%         tempCZ = myvec(Z);
        A(ic,ic,ipeak) = lsqnonneg(Z(ic,:)',tempX(ic,:)');
    end
end

%%%% update Z %%%%%
tempU = (C*A(:,:,1))'*(phi^-1)*(C*A(:,:,1));
tempW = (C*A(:,:,1))'*(phi^-1)*BCG(:,:,1);
for ipeak = 2:trial
    tempU = tempU + (C*A(:,:,ipeak))'*phi^-1*(C*A(:,:,ipeak));
    tempW = tempW + (C*A(:,:,ipeak))'*phi^-1*BCG(:,:,ipeak);
end
Z = tempU\tempW;

%%%% update phi %%%%
phi = zeros(chan);
for ipeak = 1:trial
    temp = BCG(:,:,ipeak) - C*A(:,:,ipeak)*Z;
    phi = phi+temp*temp';
end
phi = phi./double(trial*pacelen);



% calculate L
tempL = 0;
for ipeak = 1:trial
tempL = tempL + trace((phi^-1)*(BCG(:,:,ipeak)-C*A(:,:,ipeak)*Z)*(BCG(:,:,ipeak)-C*A(:,:,ipeak)*Z)');
end
L = double(trial*pacelen)*log(det(phi))+tempL;
disp(['iter',num2str(it),':  L=',num2str(L)]);



end
%%%% normalize %%%%
for ic = 1:nc
    amp1 = sqrt(var(C(:,ic)));
    C(:,ic) = C(:,ic)/amp1;
    amp2 = sqrt(var(Z(ic,:)));
    Z(ic,:) = Z(ic,:)/amp2;
    A(ic,ic,:) = amp1*amp2*A(ic,ic,:);
end

return
