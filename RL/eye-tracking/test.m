%
% jackie_test
%

clc 
clear all
close all

addpath('functions');
addpath('tetio');  

% *************************************************************************
%
% Initialization and connection to the Tobii Eye-tracker
%
% *************************************************************************
 
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

% *************************************************************************
%
% Display a stimulus 
%
% The first presentation of the stressor stimuli (or control)
%
% Participants will be instructed to fixate on a cross at the center of the
% screen with hand submerged in water. Gaze Data will be tracked and saved.
%
% *************************************************************************
close all;

%skip sync tests (currently an error)
Screen('Preference', 'SkipSyncTests', 1);

%open screen
[wPtr, rect] = Screen('OpenWindow', 0);

%write instructions for gaze
myText = 'Please keep gaze on screen.';
DrawFormattedText(wPtr,myText,'center',rect(4)/2,0);
Screen('Flip',wPtr);
WaitSecs(5);

% Drawing the cross
crossLength = 20;
crossWidth = 5;
crossColor = 0; 

drawFixationCross(wPtr,rect,crossLength,crossColor,crossWidth);
Screen('Flip',wPtr);
WaitSecs(5);

%instructions for water
myText = 'In the following part of the experiment, you are asked to immerse your hand in the water.\n\n\n\nThe experimenter will let you know when you are able to take your hand out of the water.\n\nOnly if you are not able to tolerate the cold water any more,\n\nare you allowed to take your hand out of the water before you are told to do so by the experimenter.\n\n\n\nHowever, please keep your hand in the water for as long as possible!\n\n\n\n[Press ANY KEY to continue.]';
DrawFormattedText(wPtr,myText,'center',rect(4)/3,0);
Screen('Flip',wPtr);
KbWait();
clear myText;

myText = 'Please inform experimenter you have reached this section.\n\n\n\nPlease keep looking at the screen as you put your hand into cold water.\n\n\n\n[Experimenter will start this section for you.]';
DrawFormattedText(wPtr,myText,'center',rect(4)/2.5,0);
Screen('Flip',wPtr);
WaitSecs(1);
KbWait();

%redraw cross
drawFixationCross(wPtr,rect,crossLength,crossColor,crossWidth);
Screen('Flip',wPtr);

%Wait 2 minutes or until key press 
% CHANGE TIME

hold on;



% *************************************************************************
%
% Start tracking and plot the gaze data read from the tracker.
%
% *************************************************************************

tetio_startTracking;

% leftEyeAll = [];
% rightEyeAll = [];
% timeStampAll = [];

pauseTimeInSeconds = 0.01;
%change to 2 minutes?
durationInSeconds = 5;

[leftEyeAll, rightEyeAll, timeStampAll] = DataCollect(durationInSeconds, pauseTimeInSeconds);

tetio_stopTracking; 
tetio_disconnectTracker; 
tetio_cleanUp;

%DisplayData(leftEyeAll, rightEyeAll);


% % Save gaze data vectors to file here using e.g:
 csvwrite('gazedataleft.csv', leftEyeAll);
 csvwrite('gazedataright.csv', rightEyeAll);
 csvwrite('gazedatatime.csv', timeStampAll);

 WaitSecs(5);
% clear Screen;

disp('Program finished.');

% *************************************************************************
%
% First iteration of task (need to add test trial before)
%
% *************************************************************************
Screen('Preference', 'SkipSyncTests', 1);
[wPtr, rect] = Screen('OpenWindow', 0);

%Instructions to start
myText = 'Welcome!.\n\n\n\nYou are about to start playing the real experiment.\n\n\n\nRemember to pay attention to your choices.\n\nOne of the upcoming trials will be randomly selected for payment.\n\n\n\nTo continue to the next slide, press ''1''.'
DrawFormattedText(wPtr,myText,'center',rect(4)/2.75,0);
Screen('Flip',wPtr);
%MAKE THIS CONDITIONAL TO NUMBER 1
KbWait();

myText = 'Block 1\n\n\n\nPress ''1'' to start.';
DrawFormattedText(wPtr,myText,'center',rect(4)/2.75,0);
Screen('Flip',wPtr);
%MAKE THIS CONDITIONAL TO NUMBER 1
KbWait();

%draw white dot in center - delay for 5 seconds 


%insert picture of task for 3 seconds 

%delay white dot for 5 seconds

%cue green square for 3 seconds

%yellow square for 0.5 seconds

%repeat 