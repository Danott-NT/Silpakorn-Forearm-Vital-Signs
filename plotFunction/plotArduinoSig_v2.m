function plotArduinoSig_v2(figID, arduinoSignal, arduinoFPS, timePerWindow, pulseRefArduinoRaw, pulseRefArduino, pulseRefFFT, pulseRefFFT_mag, pulseRefFFT_bpm, pulseRefPeak, pulseRefPeak_ind, picPath)
  f = figure('visible','off');
  pulseRefPeak_ind = pulseRefPeak_ind{1};
  subplot(4,1,1); plot(pulseRefArduinoRaw, "linewidth", 5); grid on; 
  title(['BPM-' num2str(pulseRefArduino)], 'FontSize', 10);
  t = 1:length(arduinoSignal);
  subplot(4,1,2); plot(arduinoSignal, 'r', t(pulseRefPeak_ind), arduinoSignal(pulseRefPeak_ind), 'or'); grid on; 
  totalPeriod = length(pulseRefPeak_ind)*length(arduinoSignal)/double(pulseRefPeak_ind(end)-pulseRefPeak_ind(1));
  title(['Peak-' num2str(pulseRefPeak{1}) ' ' num2str(length(pulseRefPeak_ind)) ' ' num2str(pulseRefPeak_ind(end)-pulseRefPeak_ind(1)) ' ' num2str(totalPeriod) ' ' num2str(arduinoFPS)]);

  d = diff(pulseRefPeak_ind);
  medIBI = median(diff(pulseRefPeak_ind)) * double(timePerWindow) / length(arduinoSignal);
  meanIBI = mean(diff(pulseRefPeak_ind)) * double(timePerWindow) / length(arduinoSignal);
  pulsePeak_med = 60/medIBI;
  pulsePeak_mean = 60/meanIBI;
  subplot(4,1,3); plot(d); title(['HRV / mean : ' num2str(pulsePeak_mean) ' med : ' num2str(pulsePeak_med) ' sd : ' num2str(std(d))]); ylabel('IBI')
  % sd = pulseRefSDArduino{1}
  subplot(4,1,4); plot(pulseRefFFT_bpm, pulseRefFFT_mag); grid on; 
  title(['FFT-' num2str(pulseRefFFT(2)) '/' num2str(pulseRefFFT(1))]); ylabel('magnitude'); xlabel('bpm');
  
  % [figID, mean([pulseRefArduino, pulseRefPeak{1}, pulseRefFFT(2)]), std([pulseRefArduino, pulseRefPeak{1}, pulseRefFFT(2)])]
  % arduinoSignal = (arduinoSignal-min(arduinoSignal))/(max(arduinoSignal)-min(arduinoSignal)); %norm to 0-1
  
  % colorPalette = {"b", "g", "r"};
  % for ch=1:size(windowedSig, 2)
  %   sig = windowedSig(:, ch);
  %   maxL = max([length(arduinoSignal), length(sig)]);
  %   arduinoSignal = interp1(1:length(arduinoSignal), arduinoSignal, linspace(1, length(arduinoSignal), maxL))+(max(sig)-min(arduinoSignal));
  %   sig = interp1(1:length(sig), sig, linspace(1, length(sig), maxL));
  %   windowedSig_time = interp1(1:length(windowedSig_time), windowedSig_time, linspace(1, length(windowedSig_time), maxL));

  %   subplot(3,4,ch*4-2); 
  %   plot(windowedSig_time, sig, colorPalette{ch});
  %   title(['CH' num2str(ch)]);
    
  %   subplot(3,4,ch*4-1); plot(specInd(:, 1), ySpec(:, ch), "b", specInd(:, 1)(peakInd_windowedSig{ch}), ySpec(:, ch)(peakInd_windowedSig{ch}), "or")
  %   title ([num2str(specInd(:, 1)(maxMagInd{ch})) '/' num2str(ySpec(:, ch)(maxMagInd{ch}))])
  % end
  % subplot(3,4,8); 
  % plot(specInd_sum, ySpec_sum, "b", specInd_sum(:, 1)(peakInd_sum_windowedSig{1}), ySpec_sum(peakInd_sum_windowedSig{1}), "or")
  % title ([num2str(specInd_sum(:, 1)(maxMagInd_sum{1})) '/' num2str(ySpec_sum(maxMagInd_sum{1}))])
  print(picPath);
  close f;

end