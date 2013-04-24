function bold = bold_ica_compsel( bold, comp )
%bold_ica ica for bold signal
% Usage:
%   bold = bold_ica( fnii, masknii, nc );
%
% Inputs:
%   bold:       bold structure from bold_ica
%   comp:       ic index to be removed
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
%   v0.1:   15-Apr-2013, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[imx, imy, imz, np] = size(bold.func);
data = double(squeeze(reshape(repmat(bold.mask,[1 1 1 np]).*bold.func, imx*imy*imz, 1, 1, []))');
act = bold.wts*bold.sph*data;
winv = pinv(bold.wts*bold.sph);
act(comp,:) = [];
winv(:,comp) = [];
func = reshape((winv*act)', imx, imy, imz, np);
bold.func = func;
bold.act = act;
bold.winv = winv;

end

