function [A,B,r] = cca(X,Y)

[n,p1] = size(X);
p2 = size(Y,2);

% Center the variables
X = X - repmat(mean(X,1), n, 1);
Y = Y - repmat(mean(Y,1), n, 1);

% Factor the inputs, and find a full rank set of columns if necessary
[Q1,T11,perm1] = qr(X,0);
% rankX = sum(abs(diag(T11)) > eps(abs(T11(1)))*max(n,p1));

[Q2,T22,perm2] = qr(Y,0);
% rankY = sum(abs(diag(T22)) > eps(abs(T22(1)))*max(n,p2));

% Compute canonical coefficients and canonical correlations.  For rankX >
% rankY, the economy-size version ignores the extra columns in L and rows
% in D. For rankX < rankY, need to ignore extra columns in M and D
% explicitly. Normalize A and B to give U and V unit variance.
% d = min(rankX,rankY);
[L,D,M] = svd(Q1' * Q2,0);
A = T11 \ L(:,:) * sqrt(n-1);
B = T22 \ M(:,:) * sqrt(n-1);
r = min(max(diag(D(:,:))', 0), 1); % remove roundoff errs



% Put coefficients back to their full size and their correct order
% A(perm1,:) = [A; zeros(p1-rankX,d)];
% B(perm2,:) = [B; zeros(p2-rankY,d)];



