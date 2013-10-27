function [A,S,z] = SNM(xs, xn, b, m)
% [C,N,T]=size(x);% C: channel number, N: Number of time points within a trial, T: Number of trials

if isempty(b)
    b = zeros(size(xn));
end

if size(xs,1)~=size(xn,1) || size(xs,1)~=size(b,1)
    error('Size not consistent');
end

% if nargin==1
%     m=C; iterN=10;
% elseif nargin==2
%     iterN=10;
% end


Ns = size(xs,2); Nn = size(xn,2); Nb = size(b,2);
ECovSig = xs*xs'/Ns;
ECovRes = xn*xn'/Nn + b*b'/Nb;
Wh=ECovRes^(-0.5);% the whitening matrix
xs_tilde=Wh*xs;
ECovSig=xs_tilde*xs_tilde'/Ns;
[W,D]=eig(ECovSig);
[Q,I]=sort(diag(D),'descend');
W=W(:,I(1:m))';
A=ECovRes^(0.5)*pinv(W);
S=pinv(A);
z=W*xs_tilde;

end