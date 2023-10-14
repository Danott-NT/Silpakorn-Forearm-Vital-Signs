function plotArduinoSig(figID, pulseRef, arduinoSignal, pulseRefArduinoRaw, pulseRefArduino, pulseRefSDArduino, pulseRefFFT, pulseRefPeak, peakInd_arduino, windowedSig, windowedSig_time, windowedSig_peakInd, specInd, ySpec, peakInd_windowedSig, maxMagInd, specInd_sum, ySpec_sum, peakInd_sum_windowedSig, maxMagInd_sum, picPath)
  f = figure('visible','off');
  colorPalette = {"b", "g", "r"};
  for ch=1:size(windowedSig, 2)
    sig = windowedSig(:, ch);
    maxL = max([length(arduinoSignal), length(sig)]);
    sig = interp1(1:length(sig), sig, linspace(1, length(sig), maxL));
    windowedSig_time = interp1(1:length(windowedSig_time), windowedSig_time, linspace(1, length(windowedSig_time), maxL));

    subplot(3,3,ch*3-2); 
    d = diff(windowedSig_peakInd{ch});
    plot(windowedSig_time, sig, colorPalette{ch});
    title(['CH' num2str(ch) ' ' num2str(length(windowedSig_peakInd{ch})) ' ' num2str(windowedSig_peakInd{ch}(end) - windowedSig_peakInd{ch}(1)) ' ' num2str(std(d))]);
    
    subplot(3,3,ch*3-1); plot(specInd(:, 1), ySpec(:, ch), "b", specInd(:, 1)(peakInd_windowedSig{ch}), ySpec(:, ch)(peakInd_windowedSig{ch}), "or")
    title ([num2str(specInd(:, 1)(maxMagInd{ch})) '/' num2str(ySpec(:, ch)(maxMagInd{ch}))])
  end
  subplot(3,3,6); 
  plot(specInd_sum, ySpec_sum, "b", specInd_sum(:, 1)(peakInd_sum_windowedSig{1}), ySpec_sum(peakInd_sum_windowedSig{1}), "or")
  title ([num2str(specInd_sum(:, 1)(maxMagInd_sum{1})) '/' num2str(ySpec_sum(maxMagInd_sum{1}))])
  print(picPath);
  close f;

  % f = figure('visible','off');
  % t = 1:length(arduinoSignal);
  % peakInd_arduino = peakInd_arduino{1};
  % subplot(3,4,1); plot(t, arduinoSignal, 'k', t(peakInd_arduino), arduinoSignal(peakInd_arduino), 'or'); grid on; set(gca,'XTick',[0:500:length(arduinoSignal)]);
  % title(['frame-' num2str(figID) ' Ref-' num2str(pulseRef)], 'FontSize', 10);
  % subplot(3,4,5); plot(pulseRefArduinoRaw, 'r', "linewidth", 5); grid on; title(['A-' num2str(pulseRefArduino) '/F-' num2str(pulseRefFFT(2)) '/P-' num2str(pulseRefPeak{1})]);
  % subplot(3,4,9); plot(diff(peakInd_arduino), '--o'); grid on; title(['HRV/SD-' num2str(pulseRefSDArduino)]); ylabel('samples'); xlabel('peak');

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
  % print(picPath);
  % close f;

end