%% load EGI raw data and perform fmrib_fastr, selection, filt and rereference
clear all; close all; clc;
eeglab;
%% input params
[filename, pathname, filterindex] = uigetfile('*.raw', 'Pick a EGI raw file');
[pathstr, name, ext] = fileparts(filename);
inName = name;
outPath = pathname;
foo = input('input [sliceNum, interFolds, arWin, ANC]: ');
if isempty(foo)
    sliceNum = 30; interFolds = 10; arWin = 30; ANC = 0;
else
    sliceNum = foo(1); interFolds = foo(2); arWin = foo(3); ANC = foo(4); 
end
foo = input('input [trSelMin, trSelMax, chSelR]: ');
if isempty(foo)
    trSel = [4 350]; chSelR = 0.5;
else
    trSel = [foo(1), foo(2)]; chSelR = foo(3);
end
foo = input('input [loCutOff, hiCutOff]: ');
if isempty(foo)
    loCutOff = 0.5; hiCutOff = 50;
else
    loCutOff = foo(1); hiCutOff = foo(2);
end
%% load EGI raw data
EEG = pop_readegi( [pathname filename] );
EEG.setname = inName;
%% add chanlocs and remove baseline
EEG.chanlocs=pop_chanedit(EEG.chanlocs, 'load',{ '/Users/hxs/Documents/Study/Research/Analysis/GSN-HydroCel-256.sfp', 'filetype', 'autodetect'});
%% TREV evaluation
[EEG,trIndex] = pop_selectevent(EEG, 'type', 'TREV');
eventTR = EEG.event(trIndex);
latency = zeros(1, length(eventTR));
for i = 1:length(eventTR)
    latency(i) = eventTR(i).latency;
end
[C,IA,IC] = unique(diff(latency));
if length(C) > 1
    disp('Error: Missing TR')
    for i = 1:length(IA)
        fprintf(1, 'Num = %d, EventNum = %d, TR Length = %d\n', i, trIndex(IA(i)), C(i));
    end
    str = input('Accept? Y/N [Y]: ', 's');
    if strcmp(str, 'N')
        return
    else
        str = input('Select TR Length Num: ', 's');
        trLength = C(str2double(str));
    end
else
    trLength = C(1);
end
%% add slice event
if rem(trLength, sliceNum)
    disp('Error: Not evenly distributed slices')
    fprintf(1, 'TR: %d\n', trLength);
    fprintf(1, 'sliceNum: %d\n', sliceNum);
    return
end
sliceInt = trLength/sliceNum;
sliceEvent1st = length(EEG.event)+1;
for i = 1:length(eventTR)
    for j = 1:sliceNum
        sliceEvent = struct('type', 'Slice',...
                            'latency', eventTR(i).latency+sliceInt*(j-1),...
                            'urevent', []);
        EEG.event = [EEG.event, sliceEvent];
    end
end
%% fmrib_fastr
[EEG, command] = pop_fmrib_fastr(EEG,[],interFolds,arWin,'Slice',...
                                 1,ANC,[],[],[],[],[],'auto');
EEG.setname = [inName '_fastr'];
%% delete slice event and save EEG
EEG.event(sliceEvent1st:end) = [];
EEG = pop_saveset(EEG, 'filename', [inName '_fastr'], 'filepath', outPath);

%% data clean up
chanR = zeros(1,EEG.nbchan);
for i = 1:EEG.nbchan
    chanR(i) = EEG.chanlocs(i).radius;
end
% select data point and channel range
EEG = pop_select(EEG, 'point', [eventTR(trSel(1)).latency eventTR(trSel(2)).latency+trLength],...
                      'channel', find(chanR<chSelR));
% remove baseline
EEG.data = rmbase( EEG.data );
% reref to average
EEG = pop_reref( EEG, [] );
% filt
[EEG, com, b] = pop_eegfiltnew(EEG, loCutOff, 0);
[EEG, com, b] = pop_eegfiltnew(EEG, 0, hiCutOff);
% save new set
EEG.setname = [inName '_fastr_sel_filted'];
EEG = pop_saveset(EEG, 'filename', [inName '_fastr_sel_filted'], 'filepath', outPath);
%% clean up
clear all; close all;