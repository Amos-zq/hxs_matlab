function checkerboard_stim_trig( shape,trialTime,trialNum,blockNum,cbSize )
%Checkerboard stimulus
% Syntax:  
%     checkerboard_stim(shape, trialTime, trialNum, blockNum, cbSize)
%
% Inputs:
%     shape     : 'circ' or 'rect'
%     trialTime : trial duration in second
%     trialNum  : number of trials
%     blockNum  : number of blocks
%     cbSize    : checkerboard size, percentage of the screen size
%
% Outputs:
%     
%
% Example:
%     
%
% Other m-files required: psychtoolbox
% Subfunctions: none
% MAT-files required: none
%
% Author: Xiaoshan Huang, xiaoshanhuang@gmail.com
%
% Versions:
%     v0.1: , orignal
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    shape = 'rect';
end

if nargin < 2
    trialTime = 1;
end
    
if nargin < 3
    trialNum = 50;
end

if nargin < 4
    blockNum = 1;
end

if nargin < 5
    cbSize = 0.8;
end

totalTime = trialTime*trialNum*blockNum

%% System config
warning('off','MATLAB:dispatcher:InexactMatch');
KbName('UnifyKeyNames');
onExit='execution halted by experimenter';

AssertOpenGL;

%% Trigger Config
triggerOut = true;

if triggerOut,
    IOPort('CloseAll');
    [comportHandle errmsg] = IOPort('OpenSerialPort', 'COM2', 'BaudRate=115200');
    IOPort('Flush', comportHandle);
end
getSensorResult = [1 3 0];
setThreshold = [3 5 0 0 0];
setEvent = [4 6 0 3 10 0];

%%

try
    myScreen = max(Screen('Screens'));
    [win,winRect] = Screen(myScreen,'OpenWindow');
    Screen('TextSize', win, 100);
    
    [width, height] = RectSize(winRect);
    Screen('FillRect',win,[0 0 0]);
    
    switch shape
        case 'rect'
            sidelength = 200;
            numCheckers =  ceil([width; height] ./ sidelength);
            % make an atomic checkerboard
            miniboard = eye(2,'uint8') .* 255;
            % repeat it in half of x,y, since it's 2x2
            checkerboard_heads = repmat(miniboard, ceil(0.5 .* numCheckers))';
            % invert for the other cycle
            checkerboard_tails = 255 - checkerboard_heads;
            % scale the images up
            checkerboard_heads = imresize(checkerboard_heads,sidelength,'box');
            checkerboard_tails = imresize(checkerboard_tails,sidelength,'box');
        case 'circ'
            % make circular checkerboard
            sigma = 12;
            spokes = 12;
            SUP = height;
            [checkerboard_heads, checkerboard_tails] = ...
                make_circular_checkerboard_pattern(sigma,spokes,SUP);
        otherwise
            disp('Shape not supported');
            return
    end
    % make textures clipped to screen size
    blankSize = (1-cbSize)/2;
    cropArea = round(height*blankSize):round(height*(1-blankSize));
    texture(1) = Screen('MakeTexture', win, checkerboard_heads(cropArea,cropArea));
    texture(2) = Screen('MakeTexture', win, checkerboard_tails(cropArea,cropArea));
    % don't need those anymore
    clear checkerboard_*; 
    
    % Define refresh rate.
    ifi = Screen('GetFlipInterval', win);
        
    %% Trigger Threshold
    Screen('DrawTexture',win,texture(1));
    Screen('Flip', win);
    WaitSecs(0.5);
    IOPort('Write', comportHandle, uint8(getSensorResult));
    msg = IOPort('Read', comportHandle, 1, 5);
    sensorWhite = msg(4)*256 + msg(5);
    Screen('DrawTexture',win,texture(2));
    Screen('Flip', win, 0, 2);
    WaitSecs(0.5);
    IOPort('Write', comportHandle, uint8(getSensorResult));
    msg = IOPort('Read', comportHandle, 1, 5);
    sensorBlack = msg(4)*256 + msg(5);
    trigThreshold = 0.8*(sensorWhite - sensorBlack)
    setThreshold(4) = fix(trigThreshold / 256);
    setThreshold(5) = mod(trigThreshold, 256);
    IOPort('Write', comportHandle, uint8(setThreshold));
    
    % Wait for key press ('s') to start
    startKey = KbName('s');
    [keyIsDown, secs, keyCode] = KbCheck;
    while ~keyCode(startKey)
        [keyIsDown, secs, keyCode] = KbCheck;
        assert(~keyCode(KbName('Escape')),onExit);
    end
    
    expStart = GetSecs;
    for i = 1:blockNum
        blockStart = GetSecs;
        %% Flash
        for j = 1:trialNum
            % Random
%             WaitSecs(0.05*(1+rand(1)));
            if triggerOut
                IOPort('Write', comportHandle, uint8(setEvent));
            end
            Screen('Flip', win, 0, 2);
            while GetSecs-blockStart<trialTime*j,
                [keyIsDown, secs, keyCode] = KbCheck;
                assert(~keyCode(KbName('Escape')),onExit);
            end
        end
    end
    
    totalTime = GetSecs - expStart
    Screen('CloseAll');
    IOPort('CloseAll');
    disp('done');
    warning('on','MATLAB:dispatcher:InexactMatch');


catch
    Screen('CloseAll');
    IOPort('CloseAll');
    psychrethrow(psychlasterror);

end

end

