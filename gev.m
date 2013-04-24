function R = gev( data, spmap )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
[Y,L] = max((data'*spmap).^2,[],2);
sigma_u = sum((sum(data.^2, 1)-sum(spmap(:,L).*data,1).^2), 2);
R = 1 - sigma_u/sum(sum(data.^2, 1), 2);
R = sqrt(R);
end

