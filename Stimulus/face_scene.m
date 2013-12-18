function face_scene( triggerOut, blocks )
%Crossmodulation frequency face and scene stimuli
% Syntax:  
%     face_scene
%
% Inputs:
%     blocks: selected blocks 
%
% Outputs:
%     
%
% Example:
%     face_scene
%
% Other m-files required: psychtoolbox
% Subfunctions: none
% MAT-files required: none
%
% Author: Xiaoshan Huang, xiaoshanhuang@gmail.com
%
% Versions:
%    v0.1: 2013-12-18 15:59, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Stimuli params
nbCat        = 2; % 1-face, 2-scene
nbCondPerCat = 3; % 1-fear/threat, 2-neutral, 3-happy/pleasant
nbCond       = nbCondPerCat^nbCat;
nbRepPerCond = 2;
nbBlock      = nbCond*nbRepPerCond;
nbTrial      = 12;

freqFace  = 6;
freqScene = 1;

timeStim = 3;
timeISI  = 2.45:0.1:3.55; % ISI 2450~3550ms, 12 steps, mean 3000ms

trigExpStart   = 200;
trigExpEnd     = 201;
trigBlockStart = 100;
trigBlockEnd   = 120;


if nargin < 1
    triggerOut = 0;
end

if nargin < 2
    blocks = 1:nbBlock;
end

%% System config
KbName('UnifyKeyNames');
onExit='execution halted by experimenter';
RestrictKeysForKbCheck([KbName('Space'),KbName('Escape')]);

AssertOpenGL;

%% Trigger Config

if triggerOut,
    config_io;
    % triggerPort = 'E800';
    triggerPort = 'C0C0';
    outp(hex2dec(triggerPort),0);
end

%% Load all stimuli
faces = {};     % {cond, gen, pics}
filepath = [pwd '/Face/fear/female/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    faces{1,1,i} = imread([filepath filename(i).name]);
end
filepath = [pwd '/Face/fear/male/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    faces{1,2,i} = imread([filepath filename(i).name]);
end
filepath = [pwd '/Face/neutral/female/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    faces{2,1,i} = imread([filepath filename(i).name]);
end
filepath = [pwd '/Face/neutral/male/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    faces{2,2,i} = imread([filepath filename(i).name]);
end
filepath = [pwd '/Face/happy/female/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    faces{3,1,i} = imread([filepath filename(i).name]);
end
filepath = [pwd '/Face/happy/male/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    faces{3,2,i} = imread([filepath filename(i).name]);
end

scenes = {};    % {cond, pics}
filepath = [pwd '/Scene/threat/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    scenes{1,i} = imread([filepath filename(i).name]);
end
filepath = [pwd '/Scene/neutral/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    scenes{2,i} = imread([filepath filename(i).name]);
end
filepath = [pwd '/Scene/pleasant/'];
filename = dir(filepath);
filename([filename.isdir]) = [];
for i = 1:length(filename)
    scenes{3,i} = imread([filepath filename(i).name]);
end

%% Generate stimuli sequence
% rng(0,'twister');
% seq = randperm(nbBlock);
% conds = zeros(nbCat+1,nbBlock);
% [condFace, condScene] = meshgrid(1:nbCondPerCat, 1:nbCondPerCat);
% 
% for rep = 1:nbRepPerCond
%     conds(1:nbCat,(rep-1)*nbCond+1:rep*nbCond) = [condFace(:) condScene(:)]';
%     conds(nbCat+1,(rep-1)*nbCond+1:rep*nbCond) = rep; 
% end
% conds = conds(:,seq);
conds = [
     2     1     3     1     3     3     2     3     2     2     1     1     2     3     2     3     1     1;   % face: 1-fear, 2-neutral, 3-happy
     3     3     1     2     1     2     2     2     2     3     1     2     1     3     1     3     1     3;   % scene: 1-threat, 2-neutral, 3-pleasant
     1     1     2     2     1     2     2     1     1     2     1     1     1     2     2     1     2     2];  % condition block index, 2 blocks for each condition
% seqs = zeros(nbTrial, nbBlock);
% for blk = 1:nbBlock
%     rng(blk,'twister');
%     seqs(:,blk) = randperm(nbTrial)';
% end

seqs = [
     3     2    11     8    11     4     8     6     1     2     7     5    11     4     3     4     9    10;
     6     7     9    12    10     5     1     5     9    10     2     8     2    10     9     9    12     4;
     5    10     7     6     3     2     9     8     4     9    12    10    12     5     2    12     4     9;
     7     9     8     9     1     9     3     7     5     7     1     1     6     8    10     8     3     2;
     4     6     3    10     9     8    10    10    11     6     5     3     7    12    12     6     1    11;
     8     5    10     2     5    11     7     9     6     5     3    11     9     7     5     1    10     1;
     9     4    12     5     8     7     6     4     8     3     6     4    10     1     8     5     2     7;
     1     1     4     4     6     6    11    12    10    11     8    12     8     9     7     2     8     6;
    11    12     1    11    12    10     4    11     7     4     4     2     1     2     4     3     7    12;
    10     3     2     1     7    12     2     3     3     8    11     7     3    11     6    11     6     5;
    12     8     5     3     2     3    12     1     2     1    10     6     4     3     1     7     5     3;
     2    11     6     7     4     1     5     2    12    12     9     9     5     6    11    10    11     8];

stimPic = {};
for blk = 1:nbBlock
    for trial = 1:nbTrial
        idx = seqs(trial, blk);
        % different face for each repetition
        rep = nbTrial/2*(conds(nbCat+1,blk)-1); % separate all faces for each repetition
        if idx > nbTrial/2  % assign half trial index as male
            stimPic{blk, 1, trial} = faces{conds(1,blk), 2, rep + idx - nbTrial/2};
        else                % assign half trial index as female
            stimPic{blk, 1, trial} = faces{conds(1,blk), 1, rep + idx};
        end
        % same scene for each repetition
        stimPic{blk, 2, trial} = scenes{conds(2,blk), idx};
    end
end

for blk = 1:nbBlock
    for trial = 1:nbTrial
        for cat = 1:nbCat
            if size(stimPic{blk, cat, trial},3) == 3
                stimPic{blk, cat, trial} = rgb2gray(stimPic{blk, cat, trial});
            end
        end
    end
end

%%
try
    myScreen = max(Screen('Screens'));
    [win,winRect] = Screen(myScreen,'OpenWindow');
    [width, height] = RectSize(winRect);
    [centerX, centerY] = WindowCenter(win);
    Priority(MaxPriority(win));
    [ifi,nrValidSamples,stddev] =Screen('GetFlipInterval', win);
    
    Screen('TextSize', win, 36);
    Screen('TextColor', win, WhiteIndex(win));
    Screen('FillRect', win, BlackIndex(win));

    HideCursor;

    fps = Screen('FrameRate', win);
    if fps == 0
        fps = round(1/ifi);
    end
    
    stimWeight = zeros(nbCat,fps);
    for frame = 1:fps
        stimWeight(1,frame) = 1/2*(1+sin(2*pi*freqFace*frame/fps));
        stimWeight(2,frame) = 1/2*(1+sin(2*pi*freqScene*frame/fps));
    end

    DrawFormattedText(win, 'Instruction, press space to continue', 'center', 'center');
    Screen('Flip', win);
    KbWait;
    WaitSecs(0.2);

    if triggerOut, outp(hex2dec(triggerPort),trigExpStart); end
    timeExpStart = GetSecs();
    if triggerOut, outp(hex2dec(triggerPort),0); end
    for blk = blocks
        % Draw texture to memory
        stimTex = {};
        step = 1; stepAll = nbTrial*nbCat;
        for trial = 1:nbTrial
            for cat = 1:nbCat
                for frame = 1:fps
                    stimTex{cat, trial, frame} = Screen('MakeTexture', win, stimWeight(cat,frame)*stimPic{blk, cat, trial});
                end
                stepPercent = [num2str(step/stepAll*100,'%3.1f') '%']; step = step + 1;
                DrawFormattedText(win, ['Loading Stimuli ' stepPercent], 'center', 'center');
                Screen('Flip', win);
            end
        end        
        % Wait
        DrawFormattedText(win, ['Block ' num2str(blk) '\n Press space when ready'], 'center', 'center');
        Screen('Flip', win);
        KbWait;
        % Block Start
        if triggerOut, outp(hex2dec(triggerPort),trigBlockStart+blk); end
        timeBlkStart = GetSecs();
        if triggerOut, outp(hex2dec(triggerPort),0); end
        for trial = 1:nbTrial
            % Prestimulus fixation
            DrawFormattedText(win, '+', 'center', 'center');
            Screen('Flip', win);
            WaitSecs(timeISI(seqs(trial, blk)));
            % Stimulus
            frame = 1;
            timeTrialStart = GetSecs();
            while(GetSecs()-timeTrialStart < timeStim)
                Screen('DrawTextures', win, [stimTex{2, trial, frame} stimTex{1, trial, frame}]);
                Screen('DrawTextures', win, [stimTex{2, trial, frame} stimTex{1, trial, frame}]);
                if triggerOut, outp(hex2dec(triggerPort),conds(1,blk)*10+conds(2,blk)); end
                vbl = Screen('Flip', win);
                frame = frame + 1; if frame > fps, frame = 1; end;
            end
            if triggerOut, outp(hex2dec(triggerPort),0); end
        end
        % Block End
        if triggerOut, outp(hex2dec(triggerPort),trigBlockEnd+blk); end
        timeBlkEnd = GetSecs();
        if triggerOut, outp(hex2dec(triggerPort),0); end
        % Clear texture memory
        for trial = 1:nbTrial
            for cat = 1:nbCat
                for frame = 1:fps
                    Screen('Close', stimTex{cat, trial, frame});
                end
            end
        end
        DrawFormattedText(win, 'Rest', 'center', 'center');
        Screen('Flip', win);
        WaitSecs(10);
    end
    if triggerOut, outp(hex2dec(triggerPort),trigExpEnd); end
    timeExpEnd = GetSecs();
    if triggerOut, outp(hex2dec(triggerPort),0); end
    timeExp = timeExpStart - timeExpEnd;
    RestrictKeysForKbCheck([]);
    ShowCursor;
    Screen('CloseAll');
    Priority(0);
    sprintf('Total experiment time: %fs',timeExp);
catch 
    RestrictKeysForKbCheck([]);
    ShowCursor;
    Screen('CloseAll');
    Priority(0);
    psychrethrow(psychlasterror);
end

