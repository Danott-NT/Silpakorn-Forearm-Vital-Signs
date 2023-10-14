function [startSample, stopSample] = splitArduinoIndexToWindow(timeSec, timeStart, timeStep, timePerWindow)
  startSec = int16(timeStart);
  stopSec = int16(timeStart+double(timePerWindow)-1);
  % if (stopSec>timeSec(end))
  %   stopSec = timeSec(end)-1;
  % endif
  startSec_last = floor(startSec);
  stopSec_last = floor(stopSec);
  if (startSec!=startSec_last || stopSec!=stopSec_last)
    startSample = find(timeSec==startSec_last)(1);
    startSample = startSample + ceil((find(timeSec==startSec_last+1)(1) - find(timeSec==startSec_last)(1))*timeStep);
    stopSample = find(timeSec==stopSec_last)(end);
    stopSample = stopSample + ceil((find(timeSec==stopSec_last+1)(1) - find(timeSec==stopSec_last)(1))*timeStep);
  else
    % startSample = find(timeSec==startSec)(1);
    % stopSample = find(timeSec==stopSec)(end);
    [m, startSample] = min(abs(timeSec-startSec));
    [m, stopSample] = min(abs(timeSec-(stopSec+1)));
    stopSample = stopSample-1;
  end
end