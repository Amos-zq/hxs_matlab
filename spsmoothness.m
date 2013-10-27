function [ P ] = spsmoothness( spmap, chanlocs, r)
%spsmoothness calculate the spatial smoothness of given spatial map,
%details see: F. Lotte and C. Guan, ?Spatially Regularized Common Spatial
%Patterns for EEG Classification? Pattern Recognition (ICPR), 2010.
% 
% Usage:
% 	P = spsmoothness( spmap, chanlocs, r);
%
% Inputs:
%     spmap    : spatial map
%     chanlocs : EEG chanlocs
%     r        : spatial smoothness radius
%
% Outputs:
%     P: spatial smoothness
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%     v0.1: 18-Mar-2013
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(spmap,1) ~= size(chanlocs,2)
    fprintf(1, 'Channel number mismatch');
    return
end

nc = size(chanlocs,2);
nmap = size(spmap,2);
chanCor = zeros(6, nc);
spmap = zscore(spmap);
for i = 1:nc
    chanCor(:,i) = [chanlocs(i).Y; chanlocs(i).X; chanlocs(i).Z; chanlocs(i).sph_theta; chanlocs(i).sph_phi; chanlocs(i).sph_radius];
end
radius = mean(chanCor(6,:),2);
G = zeros(nc,nc); D = G; K = G;
for i = 1:nc
    for j = 1:nc
%         G(i,j) = exp(-(1/2)*(sum((chanCor(1:3,i)-chanCor(1:3,j)).^2)/(r^2)));
        G(i,j) = exp(-(1/2)*(distance(chanCor(4,i),chanCor(5,i),chanCor(4,j),chanCor(5,j),radius)^2)/(r^2));
    end
end
for i = 1:nc
    D(i,i) = sum(G(i,:));
end;
K = D - G;

P = zeros(nmap,1);
for i = 1:nmap
    P(i) = spmap(:,i)'*K*spmap(:,i);
end

end