function [img, ovl] = disp_act_slc( dim, anat, act, amin, amax )
%disp_act_slc display activation map in slices 
%
% Useage:
%   [img ovl] = disp_act_slc( dim, anat, act, amin, amax )
%
% Inputs:
%   dim:    display dimension
%   anat:   anatomy image, in slice order
%   act:    activation image, in slice order
%   amin:   min activation threshold
%   amax:   max activation threshold
%
% Outputs:
%   img:    anatomy image slice by dim
%   ovl:    activation image by slice 
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   12-Mar-2013, original  
%   v0.2:   15-Apr-2013, optimize the overlay scale method
%   v0.3:   19-Apr-2013, adjust anat contrast
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[npx, npy, nslc] = size(anat);
img = zeros(npy*dim(1), npx*dim(2));
ovl = img;
if isempty(act)
    act_ovl = 0;
else
    act_ovl = 1;
    act = imresize(act, size(anat,1)/size(act,1));
    if size(anat) ~= size(act)
        fprintf('different dimension between anatomy and activation.');
        return
    end
    if isempty(amin)
        amin = [-2.5 2.5];
    end
    act( act>amin(1) & act<amin(2) ) = 0;
    if isempty(amax)
%         amax = [min(min(min(act))) max(max(max(act)))];
        amax = [-7.5 7.5];
    end
end

% roate image to display orientation, map to dim dimension
for ix = 1:size(img,2)
    for iy = 1:size(img,1)
        iax = rem(ix-1,npx)+1;
        iay = npy-rem(iy-1,npy);
        iaslc = floor((iy-1)/npy)*dim(2) + floor((ix-1)/npx)+1;
        if iaslc > nslc
            break;
        end
        img(iy,ix) = anat(iax, iay, iaslc);
        if act_ovl
            ovl(iy,ix) = act(iax, iay, iaslc);
        end
    end
end
img = mat2gray(img, [min(img(img~=0)), max(img(img~=0))]);
imshow(imadjust(img, [0.1; 0.9], [], 1.8), 'Border', 'tight', 'Colormap', gray(256));
% imshow(adapthisteq(img), 'Border', 'tight', 'Colormap', gray(256));

if nslc == 30
    imshow(imadjust(img, [0.1; 0.9], [], 0.5), 'Border', 'tight', 'Colormap', gray(256));
elseif nslc == 61
    imshow(imadjust(img, [0.1; 0.9], [], 1.8), 'Border', 'tight', 'Colormap', gray(256));
else
    imshow(adapthisteq(img), 'Border', 'tight', 'Colormap', gray(256));
end

% overlay activation map

% Number of color levels to create
nLevels = 256;

if act_ovl
    hold on
    ovlpind = gray2ind(mat2gray(-ovl.*(ovl<0), -[amin(1) amax(1)]), nLevels/2)...
            + gray2ind(mat2gray(ovl.*(ovl>0), [amin(2) amax(2)]), nLevels/2);
    ovlind = ovlpind+nLevels/2;
    ovlind(ovl<0) = nLevels/2-1 - ovlpind(ovl<0);
    ovlind = ovlind + 1;
    cmap = fireice(nLevels);
    ovlc = zeros([size(ovl),3]);
    for ix = 1:size(ovl,2)
        for iy = 1:size(ovl,1)
            ovlc(iy, ix, :) = cmap(ovlind(iy, ix), :);
        end
    end
    h = imshow(ovlc, 'Border', 'tight');
    alpha = zeros(size(ovl));
    alpha(ovl~=0) = 1;
    set(h, 'AlphaData', alpha);
    hold off
end

end
