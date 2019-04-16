%
% EyeTrackingSample
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
% Display a stimulus 
%
% For the demo this simply reads and display an image.
% Any method for generation and display of stimuli availble to Matlab could
% be inserted here, for example using Psychtoolbox or Cogent. 
%
% *************************************************************************
close all;


X = imread('TobiiDots.jpg');
img=X;

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
durationInSeconds = 1.5*1;

[leftEyeAll, rightEyeAll, timeStampAll] = DataCollect(durationInSeconds, pauseTimeInSeconds);

tetio_stopTracking; 
tetio_disconnectTracker; 
tetio_cleanUp;

DisplayData(leftEyeAll, rightEyeAll);

disp('Program finished.');
