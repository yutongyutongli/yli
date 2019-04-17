function [ExpleftEyeAll, ExprightEyeAll, ExptimeStampAll] = DataCollectExp(durationInSeconds,pauseTimeInSeconds )
%DATACOLLECT collects the data from the eye tracker
% This function is used to collect the incoming data from the tobii eye tracker.
%     
%     Input:
%         durationInSeconds: duration of the desired acquisition.
%         pauseTimeInSeconds: time lapse between readings.
%     
%     Output:
%         leftEyeAll: EyeArray corresponding to the left eye.
%         rightEyeAll:EyeArray corresponding to the right eye.
%         timeStampAll : timestamp of the readings

ExpleftEyeAll = [];
ExprightEyeAll = [];
ExptimeStampAll = [];

for i = 1:(durationInSeconds)
    
     pause(pauseTimeInSeconds);
    
    [Explefteye, Exprighteye, Exptimestamp, ExptrigSignal] = tetio_readGazeData_Exp;
    
    if isempty(Explefteye)
        continue;
    end
    
    ExpnumGazeData = size(Explefteye, 2);
    ExpleftEyeAll = vertcat(ExpleftEyeAll, Explefteye(:, 1:ExpnumGazeData));
    ExprightEyeAll = vertcat(ExprightEyeAll, Exprighteye(:, 1:ExpnumGazeData));
    ExptimeStampAll = vertcat(ExptimeStampAll, Exptimestamp(:,1));
    
end


end

