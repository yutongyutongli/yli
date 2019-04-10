function RL(observer)
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
%disp(observer);
SetCalibParams; 

disp('Starting TrackStatus');
% Display the track status window showing the participant's eyes (to position the participant).
TrackStatus; % Track status window will stay open until user key press.
disp('TrackStatus stopped');

disp('Starting Calibration workflow');
% Perform calibration
HandleCalibWorkflow(Calib);
disp('Calibration workflow stopped');

%% 0. If no argument is used, run practice and abort
if ~exist('observer', 'var')
  tetio_startTracking;
  RL_WM_Practice(randi(100));
  return;
end
%% 1. Get records from data files, if available
dataFolder = fullfile('data', num2str(observer));
if ~exist(dataFolder, 'dir')
  mkdir(dataFolder);
  fprintf('No existing data file for %d at %s. Creating.\n', observer, dataFolder);
end

[gainsBlockCount, lossBlockCount, dateBegun] = getInfoFromDataFiles(observer);
if ~isnan(dateBegun)
  fprintf('Reading existing data files from %s; experiment started on %s\n', ...
    dataFolder, dateBegun);
end
%% 2. Back up data folder (with all participant records)
if ~exist('backup', 'dir')
  mkdir('backup');
end
dateMark = datestr(datetime, 'YYYYmmDD_hhMM');
copyfile('data', fullfile('backup', sprintf('%s', dateMark)));
end