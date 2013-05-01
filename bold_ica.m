function bold = bold_ica( fnii, masknii, nc )
%bold_ica ica for bold signal
% Usage:
%   bold = bold_ica( fnii, masknii, nc );
%
% Inputs:
%   fnii:       f nii path
%   masknii:    mask nii path
%   nc:         Number of components
%
% Outputs:
%   bold.
%   func:       func
%   mask:       mask
%   wts:        ica weights
%   sph:        ica sphere
%   act:        ica activation
%   winv:       ica weight pvin
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   29-Mar-2013, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nii = load_nii(fnii);
func = double(nii.img);
[imx, imy, imz, np] = size(func);
if isempty(masknii)
    mask = ones(imx, imy, imz);
else
    nii = load_nii(masknii);
    mask = double(nii.img);
    mask = imresize(mask, size(func,1)/size(mask,1));
end
data = double(squeeze(reshape(repmat(mask,[1 1 1 np]).*func, imx*imy*imz, 1, 1, []))');
data(isnan(data)) = 0;
[wts,sph] = binica(data, 'pca', nc);
act = mapstd(wts*sph*data);
act = reshape(act', imx, imy, imz, nc);
winv = zscore(pinv(wts*sph));

bold = struct('func', func,...
              'mask', mask,...
              'wts', wts, 'sph', sph, 'act', act, 'winv', winv);
delete bias* binica* temp*

end

