function EEGOUT = loadeol(filepath, saveset)
%import *.EOL file saved by Mercury.U16 and save as EEGLAB dataset
% Syntax:  
%     EEG = loadeol(filepath, saveset)
%
% Inputs:
%     filepath: *.EOL file path
%     saveset:  save data set to same folder as input file , default = 0
%     don't save
%
% Outputs:
%     EEGOUT:   EEGLAB EEG structure
%
% Example:
%     EEG = loadeol(filepath, 1)
%
% Other m-files required: eeglab
% Subfunctions: none
% MAT-files required: none
%
% Author: Xiaoshan Huang, hxs@neuracle.cn
%
% Versions:
%     v0.1: 2013-12-03 20:49, orignal
%     v0.2: 2013-12-09 20:50, save as eeglab struct
%
% Copyright (c) 2013 Neuracle, Inc. All Rights Reserved. http://neuracle.cn/
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    [filename, pathname, filterindex] = uigetfile('*.eol; *.EOL', 'Select *.eol file');
    filepath = [pathname filename];
end

if nargin <2
    saveset = 0;    % don't save dataset
end

% Read channel number = EEG channels + trigger channel (1)
fid = fopen(filepath, 'r');
fseek(fid,50,'bof');
chan = fread(fid,1,'double');

fseek(fid,200,'bof');
data =fread(fid,[chan inf],'int32');
fclose(fid);

trigger = data(end, :);
data = data(1:end-1, :);

% Configurate triggers
trigger = [0 diff(trigger)];
trigger(trigger <= 0) = 0;
latency = find(trigger>0);
event = struct([]);
for i = 1:length(latency)
    event(i).latency = latency(i);
    event(i).type = num2str(trigger(latency(i)));
    event(i).duration = 1;
end

% import to EEG structure
[pathname, filename, EXT] = fileparts(filepath);
EEGOUT = pop_importdata(...
                    'setname', filename, ...
                    'data', data, ...
                    'dataformat', 'array', ...
                    'nbchan', chan, ...
                    'srate', 1000, ...
                    'ref', 'M1');
EEGOUT.event = event;

% save to file
if saveset
    EEGOUT = pop_saveset( EEGOUT, 'filename', EEG.setname, 'filepath', pathname);
end

end