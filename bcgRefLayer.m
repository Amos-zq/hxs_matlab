function [ EEG, regCoef, corrER ] = bcgRefLayer( EEG, BCG )
%BCG removal using reference layer method
% Syntax:  eeg  = bcgRefLayer( eeg, bcg )
%    
%
% Inputs:
%    EEG:
%    BCG: 
%
% Outputs:
%    EEG:
%
%
% Example:
%    
%
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
%
% Author: Xiaoshan Huang, xiaoshanhuang@gmail.com
%
% Versions:
%    v0.1: 2013-09-17 16:34, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize progress indicator1727_filted_EEG_reref_BRL epochs
nSteps = EEG.nbchan;
step = 0;
fprintf(1, 'bcgRefLayer: |');
strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
tic;
corrER = [];
regCoef = zeros(BCG.nbchan, EEG.nbchan);
for chan = 1:EEG.nbchan
    regCoef(:, chan) = BCG.data'\EEG.data(chan,:)';
    refSignal = (BCG.data'*regCoef(:, chan))';
    corrER(chan) = corr(EEG.data(chan,:)', refSignal');
    EEG.data(chan,:) = EEG.data(chan,:) - refSignal;
    [step, strLength] = mywaitbar(chan, EEG.nbchan, step, nSteps, strLength);
end
EEG.setname = [EEG.setname '_BRL'];
% Deinitialize progress indicator
fprintf(1, '\n');
% figure, topoplot(corrER, EEG.chanlocs);
function [step, strLength] = mywaitbar(compl, total, step, nSteps, strLength)
    progStrArray = '/-\|';
    tmp = floor(compl / total * nSteps);
    if tmp > step
        fprintf(1, [repmat('\b', 1, strLength) '%s'], repmat('=', 1, tmp - step))
        step = tmp;
        ete = ceil(toc / step * (nSteps - step));
        strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '%s %3d%%, ETE %02d:%02d'], progStrArray(mod(step - 1, 4) + 1), floor(step * 100 / nSteps), floor(ete / 60), mod(ete, 60));
    end
end

end

