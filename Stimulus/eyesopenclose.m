function eyesopenclose( trialTime,trialNum,blockNum )
%EYESOPENCLOSE eyes open and close experiment  stimulus
%
% Usage:
%     eyesopenclose( trialTime,trialNum,blockNum );
%
% Inputs:
%     trialTime : open/close duration
%     trialNum  : number of trails in a block
%     blockNum  : block number
%
% Outputs:
%	
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%	v0.1:   2013-09-03 17:50, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    trialTime = 15;
end
    
if nargin < 2
    trialNum = 1;
end

if nargin < 3
    blockNum = 10;
end

totalTime = trialTime*trialNum*(2*blockNum+1)

flipHorizontal = 0;
flipVertical = 0;

%% System config
warning('off','MATLAB:dispatcher:InexactMatch');
KbName('UnifyKeyNames');
onExit='execution halted by experimenter';

AssertOpenGL;

%% Trigger Config
triggerOut = false;

if triggerOut,
    config_io;
    triggerPort = 'E800';
    TRTrigger = 2;
    openTrigger = 1;
    closeTrigger = 3;
    outp(hex2dec(triggerPort),0);
end
%% Load audio files
audioFileNum = 16;
wavedata = {};
currentPath = pwd;
[y, freq] = audioread( [currentPath '/eyesopenclose_audio/open.wav'] );
openData = y';
nrchannels = size(openData,1);
if nrchannels < 2
    openData = [openData; openData];
    nrchannels = 2;
end
[y, freq] = audioread( [currentPath '/eyesopenclose_audio/close.wav'] );
closeData = y';
nrchannels = size(closeData,1);
if nrchannels < 2
    closeData = [closeData; closeData];
    nrchannels = 2;
end
%% Sound Init
% Perform basic initialization of the sound driver:
InitializePsychSound;

% Open the default audio device [], with default mode [] (==Only playback),
% and a required latencyclass of zero 0 == no low-latency mode, as well as
% a frequency of freq and nrchannels sound channels.
% This returns a handle to the audio device:
try
    % Try with the 'freq'uency we wanted:
    pahandle = PsychPortAudio('Open', [], [], 0, freq, nrchannels);
catch
    % Failed. Retry with default frequency as suggested by device:
    fprintf('\nCould not open device at wanted playback frequency of %i Hz. Will retry with device default frequency.\n', freq);
    fprintf('Sound may sound a bit out of tune, ...\n\n');

    psychlasterror('reset');
    pahandle = PsychPortAudio('Open', [], [], 0, [], nrchannels);
end

runMode = 1;
PsychPortAudio('RunMode', pahandle, runMode);
PsychPortAudio('FillBuffer', pahandle, closeData);

%%
try
    % Screen Init
    screens=Screen('Screens');
    screenNumber=max(screens);
    black=BlackIndex(screenNumber);
    white=WhiteIndex(screenNumber);
    [width, height]=Screen('WindowSize', screenNumber);
    wPtr = Screen('OpenWindow', screenNumber, black,[], [], 2);
    Screen('TextFont',wPtr,'Arial');
    Screen('TextSize',wPtr,round(height/10));
    % Draw Fixation
    DrawFormattedText(wPtr, '+', 'center', 'center', white, [], flipHorizontal, flipVertical);
    Screen('Flip',wPtr);
    
    % Wait for key press ('s') to start
    startKey = KbName('s');
    [keyIsDown, secs, keyCode] = KbCheck;
    while ~keyCode(startKey)
        [keyIsDown, secs, keyCode] = KbCheck;
        assert(~keyCode(KbName('Escape')),onExit);
    end
    if triggerOut, outp(hex2dec(triggerPort),TRTrigger); end
    
    expStart = GetSecs;
    % Prestimulus Fixation
    while GetSecs-expStart<trialTime*trialNum,
        [keyIsDown, secs, keyCode] = KbCheck;
        assert(~keyCode(KbName('Escape')),onExit);
    end
%     WaitSecs(1);
    if triggerOut, outp(hex2dec(triggerPort),0); end
    for i = 1:blockNum
        blockStart = GetSecs;
        for j = 1:trialNum
            % Random
%             WaitSecs(0.05*(1+rand(1)));
            % close
            if triggerOut, outp(hex2dec(triggerPort),0); end
            DrawFormattedText(wPtr, 'Close', 'center', 'center', white, [], flipHorizontal, flipVertical);
            Screen('Flip',wPtr);
            % Start audio immediately and wait for the playback to start
            PsychPortAudio('Start', pahandle, 1, 0, 1);
            if triggerOut, outp(hex2dec(triggerPort),closeTrigger); end
            % Fill buffer with next audio data
            PsychPortAudio('FillBuffer', pahandle, openData);
            while GetSecs-blockStart<trialTime*(2*j-1),
                [keyIsDown, secs, keyCode] = KbCheck;
                assert(~keyCode(KbName('Escape')),onExit);
            end
            % open
            if triggerOut, outp(hex2dec(triggerPort),0); end
            DrawFormattedText(wPtr, 'Open', 'center', 'center', white, [], flipHorizontal, flipVertical);
            Screen('Flip',wPtr);
            % Start audio immediately and wait for the playback to start
            PsychPortAudio('Start', pahandle, 1, 0, 1);
            if triggerOut, outp(hex2dec(triggerPort),openTrigger); end
            % Fill buffer with next audio data
            PsychPortAudio('FillBuffer', pahandle, closeData);
            while GetSecs-blockStart<2*trialTime*j,
                [keyIsDown, secs, keyCode] = KbCheck;
                assert(~keyCode(KbName('Escape')),onExit);
            end
            if triggerOut, outp(hex2dec(triggerPort),0); end
        end
    end
    totalTime = GetSecs - expStart
    PsychPortAudio('Stop', pahandle);
    PsychPortAudio('Close', pahandle);
    Screen('CloseAll');
    frpintf('done');
    warning('on','MATLAB:dispatcher:InexactMatch');
    
catch
    PsychPortAudio('Stop', pahandle);
    PsychPortAudio('Close', pahandle);
    Screen('CloseAll');
    psychrethrow(psychlasterror);
end

end

