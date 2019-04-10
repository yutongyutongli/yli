function results = RL_WM(observer)
clear;
clc;
% Make sure the script is running on Psychtoolbox-3:
AssertOpenGL;
% %set default values for input arguments
% if ~exist('subID','var')
%     subID=166;
% end

% %warn if duplicate sub ID
% fileName=['MLIexpSubj' num2str(subID) '.txt'];
% if exist(fileName,'file')
%     if ~IsOctave
%         resp=questdlg({['the file ' fileName 'already exists']; 'do you want to overwrite it?'},...
%             'duplicate warning','cancel','ok','ok');
%     else
%         resp=input(['the file ' fileName ' already exists. do you want to overwrite it? [Type ok for overwrite]'], 's');
%     end
%     
%     if ~strcmp(resp,'ok') %abort experiment if overwriting was not confirmed
%         disp('experiment aborted')
%         return
%     end
% end

KbName('KeyNamesWindows');
KbCheck;

%Set higher DebugLevel, so that you don't get all kinds of messages flashed
%at you each time you start the experiment:
olddebuglevel=Screen('Preference', 'VisualDebuglevel', 3);
try
    numTrial=3;
    condtable.CorrAns=nan(1,numTrial);
    condtable.CorrAns(:)=1;
    condtable.DcsDelay=nan(1,numTrial);
    condtable.DcsDelay(1:2)=3;
    condtable.DcsDelay(3)=1;
    condtable.ITIDelay=nan(1,numTrial);
    condtable.ITIDelay(1:2)=3;
    condtable.ITIDelay(3)=1;
    
    % Enable unified mode of KbName, so KbName accepts identical key names on
    % all operating systems (not absolutely necessary, but good practice):
    
    s = RandStream.create('mt19937ar','seed',sum(100*clock));
    RandStream.setGlobalStream(s);
    % whichScreen = 0;
    % thisdir = pwd;
    % Data.observer=observer;
    Data.date=datestr(now,'yyyymmddTHHMMSS');
    Data.choice=nan(1,numTrial);
    
    %disable output of keypresses to Matlab. !!!use with care!!!!!!
    %if the program gets stuck you might end up with a dead keyboard
    %if this happens, press CTRL-C to reenable keyboard handling -- it is
    %the only key still recognized.
    %ListenChar(2);
    
    %Choosing the display with the highest display number is
    %a best guess about where you want the stimulus displayed.
    %usually there will be only one screen with id = 0, unless you use a
    %multi-display setup:
    screens=Screen('Screens');
    screenNumber=max(screens);
    
    % Define black and white
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    % Open an on screen window
    [windowRect, window] = PsychImaging('OpenWindow', screenNumber, black);
    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', windowRect);
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(window);
    
    %get rid of the mouse cursor, we don't have anything to click at anyway
    HideCursor;
    
    % Preparing and displaying the welcome screen
    % We choose a text size of 24 pixels - Well readable on most screens:
    Screen('TextSize', windowRect, 36);
    
    % This is our intro text. The '\n' sequence creates a line-feed:
    myText = ['In this experiment you are asked to make choices\n' ...
        'given two rectangles on the screen.\n\n' ...
        '  Press  1  if you want to choose the left one.\n' ...
        '  Press  2  if you want to choose the right one.\n' ...
        'You will begin with ' num2str(numTrial) ' training trials\n\n' ...
        '      (Press any key to start training)\n' ];
    
    % Draw 'myText', centered in the display window:
    DrawFormattedText(windowRect, myText, 'center', 'center',white);
    
    % Show the drawn text at next display refresh cycle:
    Screen('Flip', windowRect);
    
    % Wait for key stroke. This will first make sure all keys are
    % released, then wait for a keypress and release:
    KbWait([], 3);
    % waitForBackTickt;
    % KbStrokeWait;
    
    for trial=1:numTrial
        %% Fixation oval
        Screen('FillOval', windowRect,white,[window(3)/2-20 window(4)/2-20 window(3)/2+20 window(4)/2+20]);
        Screen('Flip',windowRect);
        WaitSecs(3);
        %     Data.trialTime(trial).feedbackStartTime=datevec(now);
        %     elapsedTime=etime(datevec(now),Data.trialTime(trial).feedbackStartTime);
        
        %% Blue rect stimulus
        % Make a base Rect of 200 by 200 pixels
        baseRect = [0 0 100 200];
        % Screen X positions of our three rectangles
        squareXpos = [screenXpixels * 0.25 screenXpixels * 0.75];
        numSqaures = length(squareXpos);
        
        % Set the colors to Blue
        allColors = [0 0;0 0; 255 255];
        
        % Make our rectangle coordinates
        allRects = nan(4, 2);
        for i = 1:numSqaures
            allRects(:, i) = CenterRectOnPointd(baseRect, squareXpos(i), yCenter);
        end
        
        % Draw the rect to the screen
        Screen('FillRect', windowRect, allColors, allRects);
        
        % Flip to the screen
        [VBLTimestamp, StimulusOnsetTime, FlipTimestamp]=Screen('Flip', windowRect);
        tic;
        %these different timestamps are not exactly the same, e.g.:
        %   plot([VBLTimestamp StimulusOnsetTime FlipTimestamp tic])
        %the difference is negligible for most experiments
        
        %record response time, two methods again
        %this is just to compare between Matlab and PTB timing.
        %In your experiment, you should settle for one method -->
        %the Psychtoolbox method of using 'StimulusOnsetTime' seems to be
        %the more reliable solution, specifically on varying hardware
        %setups or under suboptimal conditions
        [resptime, keyCode] = KbWait([], 2);
        MLrt=toc;
        rt=resptime-StimulusOnsetTime;
        cc=KbName(keyCode);
        %% Decision visualization
        Data.trialTime(trial).respStartTime=datevec(now);
        elapsedTime=etime(datevec(now),Data.trialTime(trial).respStartTime);
        while elapsedTime<inf
            [keyisdown, secs, keycode, deltaSecs] = KbCheck;
            if keyisdown && (keycode(KbName('2@')) || keycode(KbName('1!')))
                elapsedTime=etime(datevec(now),Data.trialTime(trial).respStartTime);
                break
            end
            elapsedTime=etime(datevec(now),Data.trialTime(trial).respStartTime);
        end
        if keyisdown && keycode(KbName('1!'))
            Data.choice(trial)=1;
            %Data.rt(trial)=rt;
            feedbackRect = [0 0 150 250];
            fbsquareXpos = [screenXpixels * 0.25];         
        elseif keyisdown && keycode(KbName('2@'))
            Data.choice(trial)=2;
            %Data.rt(trial)=rt;
            feedbackRect = [0 0 150 250];
            % fbcolor=[255; 255; 255];
            fbsquareXpos = [screenXpixels * 0.75];
        end
        
        fbRect=CenterRectOnPointd(feedbackRect,fbsquareXpos,yCenter);
        Screen('FillRect', windowRect, [255 255 255], fbRect);
        
        baseRect = [0 0 100 200];
        allColors = [0 0; 0 0;255 255];
        % Screen X positions of our two blue rectangles
        squareXpos = [screenXpixels * 0.25 screenXpixels * 0.75];
        numSqaures = length(squareXpos);
        % Make our rectangle coordinates
        allRects = nan(4, 2);
        for i = 1:numSqaures
            allRects(:, i) = CenterRectOnPointd(baseRect, squareXpos(i), yCenter);
        end
        
        Screen('FillRect', windowRect, allColors, allRects);
        % Flip to the screen
        Screen('Flip', windowRect);
        WaitSecs(condtable.DcsDelay(trial));
        
        %% Feedback
        if Data.choice(trial)==condtable.CorrAns(trial)
            Screen('TextSize', windowRect, 70);
            RewardText = ['+ $5'];
        else
            Screen('TextSize', windowRect, 70);
            RewardText = ['+ $0'];
        end
        DrawFormattedText(windowRect, RewardText, 'center', 'center',[255 255 0]);
        Screen('Flip', windowRect);
        WaitSecs(3);
        
        %% ITI delay
        Screen('FillRect', windowRect, [0 0 0], window);
        Screen('Flip', windowRect);
        WaitSecs(condtable.ITIDelay(trial));
        % KbStrokeWait;
    end % end of trials loop
    sca;
    % Clear the screen
    
catch
    % This section is executed only in case an error happens in the
    % experiment code implemented between try and catch...
    ShowCursor;
    sca; %or sca
    ListenChar(0);
    Screen('Preference', 'VisualDebuglevel', olddebuglevel);
    %output the error message
    psychrethrow(psychlasterror);
end
end

function waitForBackTick;
while 1
    [keyisdown, secs, keycode, deltaSecs] = KbCheck;
    if keyisdown && keycode(KbName('5%'))==1
        break
    end
end
end
