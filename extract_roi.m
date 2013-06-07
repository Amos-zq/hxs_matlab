function roi_tc = extract_roi( fnii, masknii )
%extract_roi extract roi time course
% Usage:
%   roi_tc = extract_roi( fnii, masknii );
%
% Inputs:
%   fnii:       f nii path
%   masknii:    roi mask nii path
%
% Outputs:
%   roi_tc:     roi time course
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   06-Jun-2013, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nii = load_nii(fnii);
func = double(nii.img);
nii = load_nii(masknii);
mask = double(nii.img);
data = double(squeeze(reshape(func, imx*imy*imz, 1, 1, [])));
roi_tc = mean(data(mask>0,:),1);
end

