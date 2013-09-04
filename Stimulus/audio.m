function audio( seq,trialTime,trialNum,blockNum )
%AUDIO AEP experiment audio stimulus
%
% Usage:
%     audio( seq,trialTime,trialNum,blockNum );
%
% Inputs:
%     seq       : audio sequence
%     trialTime : trial duration
%     trialNum  : number of trails in a block
%     blockNum  : block number
%
% Outputs:
%	
%
% Author: Huang Xiaoshan, xiaoshanhuang@gmail.com
%
% Versions:
%	v0.1:   2013-09-03 17:42, orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    seq = 'rand';
end

if nargin < 2
    trialTime = 1;
end
    
if nargin < 3
    trialNum = 15;
end

if nargin < 4
    blockNum = 10;
end

totalTime = trialTime*trialNum*(2*blockNum+1)

%% System config
warning('off','MATLAB:dispatcher:InexactMatch');
KbName('UnifyKeyNames');
onExit='execution halted by experimenter';

AssertOpenGL;

%% Trigger Config
triggerOut = true;

if triggerOut,
    config_io;
    triggerPort = 'E800';
    TRTrigger = 2;
    stimTrigger = 1;
    outp(hex2dec(triggerPort),0);
end
% %% Load audio files
% audioFileNum = 16;
% wavedata = {};
% currentPath = pwd;
% for i = 1:audioFileNum
%     % Read WAV file from filesystem:
%     filePath = [currentPath '/Tones/TONE' num2str(i) '.WAV'];
%     [y, freq] = audioread( filePath );
%     wavedata{i} = y';
%     nrchannels = size(wavedata{i},1);
%     % Make sure we have always 2 channels stereo output.
%     if nrchannels < 2
%         wavedata{i} = [wavedata{i} ; wavedata{i}];
%         nrchannels = 2;
%     end
% end
%% Generate audio
freq = 48000;
nrchannels = 2;
duration = 0.25;
amplitude = 1;
wavedata = {};
audioFileNum = 16;
baseFreq = 100;
interFreq = 32;
for i = 1:audioFileNum
    temp = sin(linspace(0, duration*(baseFreq+interFreq*i)*2*pi, round(duration*freq)));
    wavedata{i} = [temp; temp];
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
PsychPortAudio('FillBuffer', pahandle, wavedata{round(audioFileNum/2)});

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
    DrawFormattedText(wPtr, '+', 'center', 'center', white);
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
%     while GetSecs-expStart<trialTime*trialNum,
%         [keyIsDown, secs, keyCode] = KbCheck;
%         assert(~keyCode(KbName('Escape')),onExit);
%     end
    WaitSecs(1);
    if triggerOut, outp(hex2dec(triggerPort),0); end
    for i = 1:blockNum
        blockStart = GetSecs;
        for j = 1:trialNum
            % Random
%             WaitSecs(0.05*(1+rand(1)));
            if triggerOut, outp(hex2dec(triggerPort),0); end
            % Start audio immediately and wait for the playback to start
            PsychPortAudio('Start', pahandle, 1, 0, 1);
            if triggerOut, outp(hex2dec(triggerPort),stimTrigger); end
            switch seq
                case 'rand'
                    audioIndex = randi(audioFileNum);
                case 'fix'
                    audioIndex = round(audioFileNum/2);
                otherwise
                    audioIndex = randi(audioFileNum);
            end
            % Fill buffer with next audio data
            PsychPortAudio('FillBuffer', pahandle, wavedata{audioIndex});
            while GetSecs-blockStart<trialTime*j,
                [keyIsDown, secs, keyCode] = KbCheck;
                assert(~keyCode(KbName('Escape')),onExit);
            end
            if triggerOut, outp(hex2dec(triggerPort),0); end
        end
        %% Fixation
        while GetSecs-blockStart<trialTime*trialNum*2
            [keyIsDown, secs, keyCode] = KbCheck;
            assert(~keyCode(KbName('Escape')),onExit);
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

