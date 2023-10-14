function plotRawSig(figID, pulseRef, arduinoSignal, peakInd_arduino, windowedSig, windowedSig_nor, windowedSig_time, picPath)
  f = figure('visible','off');
  % t = 1:length(arduinoSignal);
  t = linspace(windowedSig_time(1), windowedSig_time(end), length(arduinoSignal));
  peakInd_arduino = peakInd_arduino{1};
  subplot(411); plot(t, arduinoSignal, 'k', t(peakInd_arduino), arduinoSignal(peakInd_arduino), 'or'); grid on;
  xlim([int16(t(1)) int16(t(end))]); set(gca,'XTick',[int16(t(1)):5:int16(t(end))]);
  title(['frame-' num2str(figID) ' Ref-' num2str(pulseRef)], 'FontSize', 10);
  % subplot(3,4,5); plot(pulseRefArduinoRaw, 'r', "linewidth", 5); grid on; title(['A-' num2str(pulseRefArduino) '/F-' num2str(pulseRefFFTArduino) '/P-' num2str(pulseRefPeakArduino{1})]);
  % sd = pulseRefSDArduino{1}
  % subplot(3,4,9); plot(diff(peakInd_arduino), '--o'); grid on; title(['HRV/SD-' num2str(pulseRefSDArduino{1})]); ylabel('samples'); xlabel('peak');

  % arduinoSignal = (arduinoSignal-min(arduinoSignal))/(max(arduinoSignal)-min(arduinoSignal)); %norm to 0-1
  
  colorPalette = {"b", "g", "r"};
  for ch=1:size(windowedSig, 2)
    subplot(4,1,ch+1); 
    plotyy(windowedSig_time, windowedSig(:, ch), windowedSig_time, windowedSig_nor(:, ch));
    title(colorPalette{ch})
    % plot(windowedSig_time, windowedSig(:, ch), colorPalette{ch});
    % title(['CH' num2str(ch)]);
    
  %   subplot(3,4,ch*4-1); plot(specInd(:, 1), ySpec(:, ch), "b", specInd(:, 1)(peakInd_windowedSig{ch}), ySpec(:, ch)(peakInd_windowedSig{ch}), "or")
  %   title ([num2str(specInd(:, 1)(maxMagInd{ch})) '/' num2str(ySpec(:, ch)(maxMagInd{ch}))])
  end
  % subplot(3,4,8); 
  % plot(specInd_sum, ySpec_sum, "b", specInd_sum(:, 1)(peakInd_sum_windowedSig{1}), ySpec_sum(peakInd_sum_windowedSig{1}), "or")
  % title ([num2str(specInd_sum(:, 1)(maxMagInd_sum{1})) '/' num2str(ySpec_sum(maxMagInd_sum{1}))])
  print(picPath);
  close f;

end