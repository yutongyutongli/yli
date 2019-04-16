%caliberate for fixation 1 min
clc;
clear;
observer = 3;
%% Calibrating Eye Tracker
addpath('eye-tracking');
addpath('eye-tracking/functions');
addpath('eye-tracking/tetio');  
% *************************************************************************
%
% Initialization and connection to the Tobii Eye-tracker
%
% *************************************************************************
%disp(observer);
disp('Initializing tetio...');

tetio_init();
% Set to tracker ID to the product ID of the tracker you want to connect to.
trackerId = 'TX060-204-02400064';

%   FUNCTION "SEARCH FOR TRACKERS" IF NOTSET
if (strcmp(trackerId, 'NOTSET'))
	warning('tetio_matlab:EyeTracking', 'Variable trackerId has not been set.'); 
	disp('Browsing for trackers...');

	trackerinfo = tetio_getTrackers();
	for i = 1:size(trackerinfo,2)
		disp(trackerinfo(i).ProductId);
    end
	tetio_cleanUp();
	error('Error: the variable trackerId has not been set. Edit the EyeTrackingSample.m script and replace "NOTSET" with your tracker id (should be in the list above) before running this script again.');
end

fprintf('Connecting to tracker "%s"...\n', trackerId);
tetio_connectTracker(trackerId)
	
currentFrameRate = tetio_getFrameRate;
fprintf('Frame rate: %d Hz.\n', currentFrameRate);

% *************************************************************************
%
% Calibration of a participant
%
% *************************************************************************

SetCalibParams; 

disp('Starting TrackStatus');
% Display the track status window showing the participant's eyes (to position the participant).
TrackStatus; % Track status window will stay open until user key press.
disp('TrackStatus stopped');

disp('Starting Calibration workflow');
% Perform calibration
HandleCalibWorkflow(Calib);
disp('Calibration workflow stopped');

%%
thisdir = pwd;
% make a data directory if necessary
if ~isdir(fullfile(thisdir,'data'))
    disp('Making data directory');
    mkdir('Calib');
end
% make an observer directory if necessary
datadirname = fullfile(thisdir,'Calib',num2str(observer));
if ~isdir(datadirname)
    mkdir(datadirname);
end
Data.filename=fullfile(datadirname,['Calib' num2str(observer)]);

tetio_startTracking;
screens=Screen('Screens');
    screenNumber=max(screens);
    % Define colors
    white = WhiteIndex(screenNumber);
    black = BlackIndex(screenNumber);
    [windowRect, window] = PsychImaging('OpenWindow', screenNumber, black);
    % Get the size of the on screen window
    [screenXpixels, screenYpixels] = Screen('WindowSize', windowRect);
    % Get the centre coordinate of the window
    [xCenter, yCenter] = RectCenter(window);
    Screen('FillOval', windowRect,white,[window(3)/2-20 window(4)/2-20 window(3)/2+20 window(4)/2+20]);
    Screen('Flip',windowRect);
    
    Data.trialTime.FbackStartTime=datevec(now);    
    [leftEyeAll, rightEyeAll, timeStampAll] = DataCollect(10,0.01);
%     if ( etime(datevec(now),Data.trialTime(trial).FbackStartTime) < 10)
%         [lefteye, righteye, timestamp, trigSignal] = tetio_readGazeData;
%         
%         if isempty(lefteye)
%             continue;
%         end
%         
%         numGazeData = size(lefteye, 2);
%         leftEyeAll = vertcat(leftEyeAll, lefteye(:, 1:numGazeData));
%         rightEyeAll = vertcat(rightEyeAll, righteye(:, 1:numGazeData));
%         timeStampAll = vertcat(timeStampAll, timestamp(:,1));
%         
%     end
    Data.gazeL = leftEyeAll;
    %Data.Right(trial)={rightEyeAll};
    Data.gazeR = rightEyeAll;
    %Data.TimeStamp(trial)={timeStampAll};
    Data.gazeT = timeStampAll;
    
    WaitSecs(10);
    %KbWait([],3);
    Screen('CloseAll')
    save(Data.filename,'Data')
    
  tetio_stopTracking;
  tetio_disconnectTracker;
  tetio_cleanUp;