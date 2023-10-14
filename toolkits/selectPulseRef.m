function [pulseRef, checkWindowInd] = selectPulseRef(pulseRefRaw, pulseRefFFT, pulseRefPeak, lastPulseRefRaw, windowIndex)
  pulseRefPeak = pulseRefPeak{1};
  changeFlag = 0;
  checkWindowInd = 0;

  x = [pulseRefRaw, pulseRefFFT, pulseRefPeak];
  if (std(x) > 4)
    d = [abs(pulseRefRaw-pulseRefFFT), abs(pulseRefRaw-pulseRefPeak), abs(pulseRefPeak-pulseRefFFT)];
    [minD, minDindex] = min(d);
    if (minDindex == 1)
      pulseRef = mean([pulseRefRaw, pulseRefFFT]);
    elseif (minDindex == 2)
      pulseRef = mean([pulseRefRaw, pulseRefPeak]);
    elseif (minDindex == 3)
      pulseRef = mean([pulseRefPeak, pulseRefFFT]);
    end
    if (minD > 10)
      pulseRef = -1;
      % pulseRef = 75;
      disp(["Cannot find proper ref. heart rate. Use arduino heart rate as ref. heart rate. Pls check later " num2str(windowIndex)])
      disp([pulseRefRaw, pulseRefFFT, pulseRefPeak])
      disp(pulseRef)
      checkWindowInd = windowIndex;
    end
  else
    pulseRef = mean(x);
  end

  % x = [pulseRefFFT, pulseRefPeak];
  % % if (abs(pulseRefFFT - pulseRefPeak)>10)
  % %   pulseRef = -1;
  % %   checkWindowInd = windowIndex;
  % % else 
  % if (abs(pulseRefFFT - pulseRefPeak)>10)
  %   pulseRef = pulseRefFFT;
  %   % pulseRef = pulseRefPeak;
  %   % pulseRef = 84;
  % else
  %   x = [pulseRefFFT, pulseRefPeak];
  %   pulseRef = mean(x);
  % end
  
end