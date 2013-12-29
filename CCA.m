function [Wx, Wy, r] = CCA( X, Y )
%CCA Canonical correlation analysis.
%
% Syntax:  
%     [Wx, Wy, r] = CCA( X, Y )
%
% Inputs:
%     X : n-by-p1 data matrix, n observations, p1 variable
%     Y : n-by-p2 data matrix, n observations, p2 variables
%
% Outputs:
%     Wx : columns of cononical coefficients for X
%     Wy : columns of cononical coefficients for Y
%     r  : correlation coefficients
%
% Example:
%     Wx, Wy, r] = CCA(X, Y);
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Xiaoshan Huang, xiaoshanhuang@gmail.com
%
% Versions:
%     v0.1: 2013-12-29 18:31, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 2
    error(message('stats:canoncorr:TooFewInputs'));
end

[n,p1] = size(X);
if size(Y,1) ~= n
    error(message('stats:canoncorr:InputSizeMismatch'));
elseif n == 1
    error(message('stats:canoncorr:NotEnoughData'));
end
p2 = size(Y,2);

% Center the variables
X = X - repmat(mean(X,1), n, 1);
Y = Y - repmat(mean(Y,1), n, 1);

C = cov([X Y]);
A = [zeros(p1,p1) ones(p1,p2);
     ones(p2,p1) zeros(p2,p2)];
B = [ones(p1,p1) zeros(p1,p2);
     zeros(p2,p1) ones(p2,p2)];

% Slove the generalized eigenvalue problem
[V,D] = eig(A.*C,B.*C);

[r,I] = sort(diag(D),'descend');
r = min(max(r,0),1);
V = V(:,I);
V = V(:,r>0);
r(r==0) = [];

Wx = V(1:p1,:);
Wy = V(p1+1:p1+p2,:);

end

