function K = sps_kernal( chanlocs, ss_r )
%EEG channel spatial smoothness kernal
%F. Lotte and C. Guan, ?Spatially Regularized Common Spatial Patterns for 
%EEG Classification,? Pattern Recognition (ICPR), 2010.
%
% Usage:
%     K = sps_kernal( chanlocs, ss_r )
%
% Inputs:
%     chanlocs:   EEG chanlocs
%     ss_r:       spatial smooth channel range
%
% Outputs:
%     K:  spatial smoothness kernal
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   01-Jul-2013, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            G(i,j) = exp(-(1/2)*(distance(chanCor(4,i),chanCor(5,i),chanCor(4,j),chanCor(5,j),radius)^2)/(ss_r^2));
        end
    end
    for i = 1:nc
        Diag(i,i) = sum(G(i,:));
    end;
    K = Diag - G;    
end

