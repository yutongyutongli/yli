function RL_WM_Practice(observer)
% Make sure the script is running on Psychtoolbox-3:
AssertOpenGL;
KbName('KeyNamesWindows');
KbCheck;
thisdir = pwd;

%Set higher DebugLevel, so that you don't get all kinds of messages flashed
%at you each time you start the experiment:
olddebuglevel=Screen('Preference', 'VisualDebuglevel', 3);
% make a data directory if necessary
if ~isdir(fullfile(thisdir,'data'))
    disp('Making data directory');
    mkdir('data');
end
% make an observer directory if necessary
datadirname = fullfile(thisdir,'data',num2str(observer));
if ~isdir(datadirname)
    mkdir(datadirname);
end
Data.filename=fullfile(datadirname,['Subj' num2str(observer)]);
filename2=fullfile(datadirname,['ValidTrials' num2str(observer)]);
figfilename= fullfile(datadirname,['Fig' num2str(observer)]);
try
    %% Data struct setup
    Data.stimulus.numTrial = 5; %enter number of valid trials expected to be accomplished
    Data.stimulus.maxExpDur = 5; %enter number of minutes allowed for experiment
    
    Data.stimulus.responseWindowDur=.5; 
    Data.stimulus.choiceDur=.5;
    Data.stimulus.ITIlong = 5;
    Data.stimulus.ITIshort = 1;
    Data.stimulus.feedbackDur=.5;
    Data.stimulus.CFIlong = 5;
    Data.stimulus.CFIshort = 1;
    Data.date=datestr(now,'yyyymmddTHHMMSS');
    Data.choice = nan;
    Data.reward = nan;
    %% Initiate data matrix; set alpha level and depth
    %data = [];
    data(1,1:2) = nan;
    alpha = 0.05; % threshold level for statistical test
    maxdepth = 4; %number from 1-4 with maximum length of recent choice sequences to test
    %% Display setup
    %Disable output of keypresses to Matlab.
    %press CTRL-C to reenable keyboard handling
    ListenChar(2);
    %get rid of the mouse cursor, we don't have anything to click at anyway
    HideCursor;
    %Choosing the display with the highest display number:
    screens=Screen('Screens');
    screenNumber=max(screens);
    % Define colors
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    % Define color for stimulus
    RectColor = [0 0; 204 204; 204 204]; % light blue
    FBackcolor = [255 255 255]; %white
    % Open an on screen window
    [windowRect, window] = PsychImaging('OpenWindow', screenNumber, black);
    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', windowRect);
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(window);
    % Make stimulus of 100 by 200 pixels size
    RectSize = [0 0 100 200];
    %% Introduction screen
    Screen('TextSize', windowRect, 36);
    % This is our intro text. The '\n' sequence creates a line-feed:
    IntroText = ['In this experiment you are asked to make choices\n' ...
        'given two rectangles on the screen.\n\n' ...
        '  Press  1  if you want to choose the left one.\n' ...
        '  Press  2  if you want to choose the right one.\n' ...
        'You will begin with ' num2str(Data.stimulus.numTrial) ' training trials\n\n' ...
        '      (Press any key to start training)\n' ];
    % Draw 'myText', centered in the display window:
    DrawFormattedText(windowRect, IntroText, 'center', 'center', white);
    % Show the drawn text at next display refresh cycle:
    Screen('Flip', windowRect);
    % Wait for key stroke. This will first make sure all keys are
    % released, then wait for a keypress and release:
    KbWait([], 3);
    Data.ExpStartTime=datevec(now);
    n = 1; % Number(n) increment, for the full choice history
    validtrials = 0; % Number of valid trials
    while etime(datevec(now),Data.ExpStartTime) < Data.stimulus.maxExpDur*60  
        if validtrials < Data.stimulus.numTrial
            trial = n;
            %% Generate correct ans with matching pannies algorithm
            [computerChoice,pComputerRight,biasInfo]=matching_pennies(data,maxdepth,alpha);
            Data.CompChoice(trial) = computerChoice;
            Data.pCompRight(trial) = pComputerRight;
            Data.biasInfo(trial,1:4) = biasInfo;
            %% Fixation oval
            Screen('FillOval', windowRect,white,[window(3)/2-20 window(4)/2-20 window(3)/2+20 window(4)/2+20]);
            Screen('Flip',windowRect);
            
            %Clock
            Data.trialTime(trial).TrialStartTime=datevec(now);
            Data.trialTime(trial).FixStartTime=datevec(now);
            
            WaitSecs(0.5);
            %Clock
            Data.trialTime(trial).FixEndTime=datevec(now);
            Data.trialTime(trial).FixDuration=etime(datevec(now),Data.trialTime(trial).FixStartTime);
            %% Scample - Blue rectangle stimulus
            % Screen X positions of our three rectangles
            squareXpos = [screenXpixels * 0.25 screenXpixels * 0.75];
            numSqaures = length(squareXpos);
            % Make our rectangle coordinates
            Rects = nan(4, 2);
            for i = 1:numSqaures
                Rects(:, i) = CenterRectOnPointd(RectSize, squareXpos(i), yCenter);
            end
            % Draw blue rectangles to the screen
            Screen('FillRect', windowRect, RectColor, Rects);
            Screen('Flip', windowRect);
            % Clock
            Data.trialTime(trial).StimStartTime=datevec(now);
            
            elapsedTime=etime(datevec(now),Data.trialTime(trial).StimStartTime);
            while elapsedTime<=Data.stimulus.responseWindowDur
                [keyisdown, secs, keycode, deltaSecs] = KbCheck;
                if keyisdown && (keycode(KbName('2@')) || keycode(KbName('1!')))
                    elapsedTime=etime(datevec(now),Data.trialTime(trial).StimStartTime);
                    break
                end
                elapsedTime=etime(datevec(now),Data.trialTime(trial).StimStartTime);
            end
            
            % Clock
            Data.trialTime(trial).StimEndTime = datevec(now);
            Data.trialTime(trial).StimDuration = elapsedTime;
            %% Decision response
            if keyisdown && keycode(KbName('1!'))
                Data.choice(trial)= 0;
                ChoiceRectSize = [0 0 150 250];
                fbsquareXpos = [screenXpixels * 0.25];
            elseif keyisdown && keycode(KbName('2@'))
                Data.choice(trial)= 1;
                ChoiceRectSize = [0 0 150 250];
                fbsquareXpos = [screenXpixels * 0.75];
            else
                Data.choice(trial)= nan;
                ChoiceRectSize = [0 0 0 0];
                fbsquareXpos = [screenXpixels * 0.75];
            end
            
            fbRect=CenterRectOnPointd(ChoiceRectSize,fbsquareXpos,yCenter);
            Screen('FillRect', windowRect, FBackcolor, fbRect);
            
            % Screen X positions of our two blue rectangles
            squareXpos = [screenXpixels * 0.25 screenXpixels * 0.75];
            numSqaures = length(squareXpos);
            % Make our rectangle coordinates
            Rects = nan(4, 2);
            for i = 1:numSqaures
                Rects(:, i) = CenterRectOnPointd(RectSize, squareXpos(i), yCenter);
            end
            Screen('FillRect', windowRect, RectColor, Rects);
            % Flip to the screen
            Screen('Flip', windowRect);
            % Clock
            Data.trialTime(trial).RespStartTime=datevec(now);
            
            WaitSecs(0.5);
            
            % Clock
            Data.trialTime(trial).RespEndTime=datevec(now);
            %% CFI
            Screen('FillOval', windowRect,white,[window(3)/2-20 window(4)/2-20 window(3)/2+20 window(4)/2+20]);
            Screen('Flip',windowRect);
            Data.trialTime(trial).CFIStartTime=datevec(now);
            RandCFI = randi(100);
            if RandCFI > 50
                WaitSecs(Data.stimulus.CFIlong);
            else
                WaitSecs(Data.stimulus.CFIshort);
            end
            Data.trialTime(trial).CFIEndTime=datevec(now);
            Data.trialTime(trial).CFIDuration=etime(datevec(now),Data.trialTime(trial).CFIStartTime);
            %% Feedback
            if isnan(Data.choice(trial))
                Screen('TextSize', windowRect, 70);
                RewardText = ['+ 0'];
                Data.reward(trial) = nan;
            else
                if Data.choice(trial)==Data.CompChoice(trial)
                    Screen('TextSize', windowRect, 70);
                    RewardText = ['+ 5'];
                    Data.reward(trial) = 1;
                else
                    Screen('TextSize', windowRect, 70);
                    RewardText = ['+ 0'];
                    Data.reward(trial) = 0;
                end
            end
            
            DrawFormattedText(windowRect, RewardText, 'center', 'center',[255 255 0]);
            Screen('Flip', windowRect);
            %Clock
            Data.trialTime(trial).FbackStartTime=datevec(now);
            
            WaitSecs(Data.stimulus.feedbackDur);
            
            %Clock
            Data.trialTime(trial).FbackEndTime=datevec(now);
            Data.trialTime(trial).FbackDuration=etime(datevec(now),Data.trialTime(trial).FbackStartTime);
            %% ITI delay
            Screen('FillOval', windowRect,white,[window(3)/2-20 window(4)/2-20 window(3)/2+20 window(4)/2+20]);
            Screen('Flip', windowRect);
            Data.trialTime(trial).ITIStartTime=datevec(now);
            RandITI = randi(100);
            if RandITI > 50
                WaitSecs(Data.stimulus.ITIlong);
            else
                WaitSecs(Data.stimulus.ITIshort);
            end
            
            %% Clock
            Data.trialTime(trial).TrialEndTime=datevec(now);
            Data.trialTime(trial).TrialDuration=etime(datevec(now),Data.trialTime(trial).TrialStartTime);
            Data.trialTime(trial).elapsedTimeExp=etime(datevec(now),Data.ExpStartTime);
            %% Add to choice and reward history
            data(1:n,1)=Data.choice(:);
            data(1:n,2)=Data.reward(:);
            %% remove invalid trials and count number of valid trials
            data=rmmissing(data);
            validtrials = sum(~isnan(data(:,1)));
            %% Data visualization
            Data.sum(n).miss=sum(isnan(Data.choice(:)));
            Data.sum(n).left=sum(Data.choice(:)==0);
            Data.sum(n).right=sum(Data.choice(:)==1);
            Data.sum(n).reward=sum(rmmissing(Data.reward(:)));      
            %% Increment, for the full choice history
            n = n+1;  
        else
            break
        end % end of trials loop
    end % end of experimen time control while loop
    Data.ExpEndTime=datevec(now);
    Data.ExpDur=etime(datevec(now),Data.ExpStartTime);
    %% plot cumulative choice and reward history graph
    Time = struct2dataset(Data.trialTime(:));
    Result = struct2dataset(Data.sum(:));
    x=Time(:,19);
    xMax=Data.ExpDur;
    reward=Result(:,4);
    Fig.One = figure;
    plot(x,reward,'color',[0.8500 0.3250 0.0980])
    hold on
    plot(x,Result(:,1),'color',[0 0.4470 0.7410])
    plot(x,Result(:,2),'color',[0.4660 0.6740 0.1880])
    plot(x,Result(:,3),'color',[0.6350 0.0780 0.1840])
    drawnow
    title('Cumulative Choice and Reward vs. Elapsed Time')
    legend('Reward','Missed','Left','Right')
    hold off
    %% Save data and graph
    save(Data.filename,'Data')
    save(filename2,'data')
    save(figfilename,'Fig')
    Screen('TextSize', windowRect, 36);
    % This is our intro text. The '\n' sequence creates a line-feed:
    IntroText = ['This is the end.\n\nYou won ' num2str(sum(data(:,2))*5) ' points.\n\nThanks for participating!'];
    % Draw 'myText', centered in the display window:
    DrawFormattedText(windowRect, IntroText, 'center', 'center', [255 255 255]);
    % Show the drawn text at next display refresh cycle:
    Screen('Flip', windowRect);
    %WaitSecs(3);
    KbWait([],3);
    %eval(sprintf('save %s.mat Data',Data.filename))
    Screen('CloseAll')
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     Below starts the algorithm function    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [computerChoice,pComputerRight,biasInfo]=matching_pennies(data,maxdepth,alpha)

% matching pennies algorithms
%    check previous choice and reward history for any choice bias
%    if no bias, pComputerRight = 0.5
%    if bias detected, picks to minimize the subject's reward
%       (pComputerRight = 1 - p(subject choosing right)
%
% inputs:
%   data: matrix with choices/outcomes from all completed trials so far
%     (:,1) - choices from session (0 Left / 1 Right)
%     (:,2) - outcomes from session (0 Unrewarded / 1 Rewarded)
%   max depth - number from 1-4 with maximum length of recent choice
%       sequences to test
%   alpha - significance threshold for two-sided tests
%
% outputs:
%   computerChoice: 0 - left / 1 - right
%   pComputerRight: probability of computer picking right on a given trial
%       (used to draw the computer's choice)
%   biasInfo: vector containing information about the computer's bias, if any
%       (1): bias detected
%           0 - no bias / 1 if left bias / 2 if right bias
%       (2): depth of bias detected
%           -1 - no bias / 0-4 - bias depth with maximum deviation
%              deviation used to bias computer choice
%       (3): magnitude of bias (pRight - 0.5)
%          0 - no bias / <0 - left bias / >0 - right bias
%       (4): which algorithm was used to bias the computer's choice
%          0 - no bias / 1 - algorithm 1 / 2 - algorithm 2

%prepping for two-sided statistical tests
testAlpha=alpha/2;

%defaults for no bias detected
pComputerRight=0.5; bias=0; biasDepth=-1; maxDeviation=0; whichAlg=0;

if isempty(data)
    biasInfo=[];
    if rand<0.5, computerChoice=1; else, computerChoice=0; end
    return;
end

choice=data(:,1)+1; %recode as 1/2
reward=data(:,2)+1;
choice(end+1)=NaN; %placeholder for current trial
reward(end+1)=NaN;

%% algorithm 1
% checks previous choice history for bias

for depth=0:maxdepth
    
    %skip testing if there aren't enough previous trials
    if length(data)<depth+1
        continue;
    end
    
    if depth==0  %count all right choices to check overall side bias
        
        countRight=sum(data(:,1)); %total number of rightward choices (X)
        countN=length(data(:,1)); %total number of trials (N)
        
    else         %look for bias on trials where recent choice history was the same
        
        %code recent choice hist for all trials
        histseq=nan(size(choice,1),1);
        for trial=depth+1:length(choice)
            seq=0;
            for currdepth=1:depth
                seq=seq+choice(trial-currdepth)*10^(currdepth-1);
            end
            histseq(trial)=seq;
        end
        
        %find all previous trials with same preceding choice sequence as current trial
        idx=find(histseq==histseq(end));
        
        %remove current trial
        idx=idx(1:end-1);
        
        %skip the statistical tests if there are no other instances of this sequence
        if isempty(idx), continue; end
        
        %count rightward choices made on trials following this sequence
        %using choices coded as 1/0
        countRight=sum(data(idx,1));
        
        %count number of previous instances of the sequence
        countN=length(idx);
        
    end
    
    %calculate probabilities using binomial CDF
    %eg. probability of choosing right X times out of N trials,
    %given the previous choice sequence
    %if this is less than alpha, we say the subject has a (conditional) side bias
    pRightBias=1-binocdf(countRight-1,countN,0.5); %p(X>=x)
    pLeftBias=binocdf(countRight,countN,0.5); %p(X<=x)
    pDeviation=(countRight/countN) - 0.5; %p(Right)
    
    %if significant bias detected, compare the deviation to current maximum
    %   since the subject may show bias for multiple depths, we use the
    %   maximum deviation from 0.5 to determine which bias the
    %   computer will exploit when making its choice
    if pRightBias<testAlpha || pLeftBias<testAlpha
        if abs(pDeviation)>abs(maxDeviation)
            maxDeviation=pDeviation;
            if maxDeviation<0, bias=1; else, bias=2; end
            biasDepth=depth;
            pComputerRight=1-(maxDeviation+0.5);
            whichAlg=1;
        end
    end
    
end

%% algorithm 2
% checks choice and reward history for bias
% this is the same as above except there is no depth 0 and the reward
% history is also included
for depth=1:maxdepth
    
    if length(data)<depth+1
        continue;
    end
    
    chistseq=nan(size(choice,1),1);
    rhistseq=nan(size(reward,1),1);
    for trial=depth+1:length(choice)
        cseq=0; rseq=0;
        for currdepth=1:depth
            cseq=cseq+choice(trial-currdepth)*10^(currdepth-1);
            rseq=rseq+reward(trial-currdepth)*10^(currdepth-1);
        end
        chistseq(trial)=cseq;
        rhistseq(trial)=rseq;
    end
    
    idx=find(chistseq==chistseq(end) & rhistseq==rhistseq(end));
    idx=idx(1:end-1);
    if isempty(idx), continue; end
    countRight=sum(data(idx,1));
    countN=length(idx);
    
    pRightBias=1-binocdf(countRight-1,countN,0.5);
    pLeftBias=binocdf(countRight,countN,0.5);
    pDeviation=(countRight/countN) - 0.5;
    
    if pRightBias<testAlpha || pLeftBias<testAlpha
        if abs(pDeviation)>abs(maxDeviation)
            maxDeviation=pDeviation;
            if maxDeviation<0, bias=1; else, bias=2; end
            biasDepth=depth;
            pComputerRight=1-(maxDeviation+0.5);
            whichAlg=2;
        end
    end
end

biasInfo=[bias biasDepth maxDeviation whichAlg];

%% make computer's choice

if rand<pComputerRight, computerChoice=1; else, computerChoice=0; end
end