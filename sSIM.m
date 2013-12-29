function [A,S,z]=sSIM(x,lambda,m,iterN)
% Matlab code for the SIgnal-to-noise ratio Maximizer (SIM) algorithm.
% SIM is an algorithm for optimizing spatial filters by maximizing the SNR
% of ERPs.
% Inputs: x: EEG signal organized by channel/time/trial; 
%         m: number of spatial filters to be optimized (default: equal to the number of channels)
%         iterN: maximum number of iterations (default: 10)
% Outputs: A: spatial patterns, each column corresponding to the spatial
%             pattern of a source
%          S: spatial filters, each row containing a spatial filter
%          z: source signals, each row containing the time course of a
%             source signal.
%
% written by Wei Wu.  Email: weiwu@neurostat.mit.edu
% Last updated on 20 July, 2011

[C,N,T]=size(x);% C: channel number, N: Number of time points within a trial, T: Number of trials
if nargin==2
    m=C; iterN=10;
elseif nargin==3
    iterN=10;
end
xs=squeeze(mean(x,3));
xn=x(:,:)-repmat(xs,1,T);
ECovSig=xs*xs'/N; % signal covariance 
% ECovSig = GraphicalLasso(xs', lambda);
ECovRes=xn*xn'/(N*T); % noise covariance
% ECovRes = GraphicalLasso(xn', lambda);
for n=1:iterN
    Wh=ECovRes^(-0.5);% the whitening matrix
    xs_tilde=Wh*xs;
%     ECovSig=xs_tilde*xs_tilde'/N;
    ECovSig = GraphicalLasso(xs_tilde', lambda);
    [W,D]=eig(ECovSig);
    [Q,I]=sort(diag(D),'descend');
    W=W(:,I(1:m))';
    A=ECovRes^(0.5)*pinv(W);
    S=pinv(A);
    z=W*xs_tilde;
%     loglike(n)=-N/2*(T*N*log(2*pi)+T*log(det(ECovRes))+1/N*trace(inv(ECovRes)*xn*xn'));% log-likelihood
    xn=x(:,:)-repmat(A*z,1,T);
    ECovRes=xn*xn'/(N*T);
%     ECovRes = GraphicalLasso(xn', lambda);
end


