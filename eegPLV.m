function plv = eegPLV( EEG, filtRange )
%Computes the Phase Locking Value (PLV) for an EEG dataset.
% Usage:
%   plv = eegPLV( EEG, filtRange )
%
% Inputs:
%   EEG         : EEGLAB data structure
%   filtRange   : filter range, e.g. [8 13]
%
% Outputs:
%   plv: Phase locking value, numTimePoints x numChannels x numChannels
%
% Reference: 
% J. P. Lachaux, E. Rodriguez, J. Martinerie, and F. J. Varela, 
% ?Measuring phase synchrony in brain signals.,? 
% Hum Brain Mapp, vol. 8, no. 4, pp. 194?208, 1999.
%
% Code Reference: pn_eegPLV.m
% http://praneethnamburi.wordpress.com/2011/08/10/plv/
% Praneeth Namburi
% Cognitive Neuroscience Lab, DUKE-NUS
% 01 Dec 2009
% Present address: Neuroscience Graduate Program, MIT
% email:           praneeth@mit.edu
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%   v0.1:   2013-07-25 17:08, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numChannels = EEG.nbchan;
numTrials = EEG.trials;
dataSelectArr = true(numTrials, 1);
numConditions = size(dataSelectArr, 2);

% disp('Filtering data...');
EEG = pop_eegfiltnew(EEG, filtRange(1), filtRange(2));
filteredData = EEG.data;

% disp(['Calculating PLV for ' mat2str(sum(dataSelectArr, 1)) ' trials...']);
for channelCount = 1:numChannels
    filteredData(channelCount, :, :) = angle(hilbert(squeeze(filteredData(channelCount, :, :))));
end
plv = zeros(size(filteredData, 2), numChannels, numChannels, numConditions);

% Initialize progress indicator
nSteps = 20;
step = 0;
fprintf(1, 'Calculating PLV: |');
strLength = fprintf(1, [repmat(' ', 1, nSteps - step) '|   0%%']);
tic

for channelCount = 1:numChannels-1
    channelData = squeeze(filteredData(channelCount, :, :));
    for compareChannelCount = channelCount+1:numChannels
        compareChannelData = squeeze(filteredData(compareChannelCount, :, :));
        for conditionCount = 1:numConditions
            plv(:, channelCount, compareChannelCount, conditionCount) = abs(sum(exp(1i*(channelData(:, dataSelectArr(:, conditionCount)) - compareChannelData(:, dataSelectArr(:, conditionCount)))), 2))/sum(dataSelectArr(:, conditionCount));
        end
    end
    for compareChannelCount = 1:channelCount
        for conditionCount = 1:numConditions
            plv(:, channelCount, compareChannelCount, conditionCount) = plv(:, compareChannelCount, channelCount, conditionCount);
        end
    end
    [step, strLength] = mywaitbar(channelCount, numChannels-1, step, nSteps, strLength);
end

% Deinitialize progress indicator
fprintf(1, '\n');

plv = squeeze(plv);



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

