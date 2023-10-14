function plotICASig(arduinoSignal, sig_ica, timePerWindow, sig_ica_peakInd, specInd_ica, ySpec_ica, peakInd_ica, magMaxInd_ica, pulseMaxInd_ica, specInd_ica_sum, ySpec_ica_sum, peakInd_ica_sum, magMaxInd_ica_sum, pulseMaxInd_ica_sum, picPath)
  % , specInd_ica_norm, ySpec_ica_norm, peakInd_ica_norm, magMaxInd_ica_norm, pulseMaxInd_ica_norm, specInd_ica_snr, ySpec_ica_snr, peakInd_ica_snr, magMaxInd_ica_snr, pulseMaxInd_ica_snr
  f = figure('visible','off');
  fontSize = 10;
  timeInd = 1:length(sig_ica(:, 1));
  ch = size(sig_ica, 2);
  arduinoSignal = (arduinoSignal-min(arduinoSignal))/(max(arduinoSignal)-min(arduinoSignal)); %norm to 0-1
  
  colorPalette = {"b", "g", "r"};
  for (i=1:ch)
    sig = sig_ica(:, i);
    timeInd = 1:length(sig);
    % maxL = max([length(arduinoSignal), length(sig)]);
    % arduinoSignal = interp1(1:length(arduinoSignal), arduinoSignal, linspace(1, length(arduinoSignal), maxL))+(max(sig)-min(arduinoSignal));
    % sig = interp1(1:length(sig), sig, linspace(1, length(sig), maxL));
    % timeInd = interp1(1:length(timeInd), timeInd, linspace(1, length(timeInd), maxL));

    subplot(ch,3,3*i-2); 
    % plot(timeInd, sig, colorPalette{i}, timeInd(sig_ica_peakInd{i}), sig(sig_ica_peakInd{i}), 'or');
    plot(timeInd, sig, colorPalette{i}, timeInd(sig_ica_peakInd{i}), sig(sig_ica_peakInd{i}), ".k", "markersize", 7);
    d = diff(sig_ica_peakInd{i});
    meanIBI = double(median(d) * timePerWindow) / length(sig);
    pulsePeak = 60/meanIBI;
    % title(['ICA Com' num2str(i) ' ' num2str(length(sig_ica_peakInd{i})) ' ' num2str(pulsePeak) ' ' num2str(std(d))], 'FontSize', 8); 
    title(['ICA Com' num2str(i) ' SD:' num2str(round(std(d)*10)/10)], 'FontSize', 8); 
    set(gca,'XTick', floor(linspace(timeInd(1), timeInd(end), 5)));
    set(gca,'YTick', round(linspace(min(sig), max(sig), 5)*10)/10);
    [~, maxPeakInd] = max(ySpec_ica(:, i)(peakInd_ica{i}));
    subplot(ch,3,3*i-1); plot(specInd_ica(:, 1), ySpec_ica(:, i), specInd_ica(:, 1)(peakInd_ica{i}(maxPeakInd)), ySpec_ica(:, i)(peakInd_ica{i}(maxPeakInd)), "or");
    title([num2str(round(pulseMaxInd_ica{i}*10)/10) ' BPM/mag ' num2str(round(magMaxInd_ica{i}*10)/10)], 'FontSize', 8); 
    set(gca,'XTick', floor(linspace(specInd_ica(1), specInd_ica(end), 5)));
    set(gca,'YTick', floor(linspace(min(ySpec_ica(:, i)), max(ySpec_ica(:, i)), 5)));
    
  end

  % sum of fft from all comp.
  subplot(336);
  [~, maxPeakInd] = max(ySpec_ica_sum(peakInd_ica_sum{1}));
  plot(specInd_ica_sum, ySpec_ica_sum, specInd_ica_sum(peakInd_ica_sum{1}(maxPeakInd)), ySpec_ica_sum(peakInd_ica_sum{1}(maxPeakInd)), "or");
  title([num2str(round(pulseMaxInd_ica_sum{1}*10)/10) ' BPM/mag ' num2str(round(magMaxInd_ica_sum{1}*10)/10)], 'FontSize', 8); 
  set(gca,'XTick', floor(linspace(specInd_ica_sum(1), specInd_ica_sum(end), 5)));
  set(gca,'YTick', floor(linspace(min(ySpec_ica_sum), max(ySpec_ica_sum), 5)));

  % sumYspec = ySpec_ica(:, 1) + ySpec_ica(:, 2) + ySpec_ica(:, 3);
  % [peakMag, peakInd, pulseRefPeak, HRVpeak] = findpeaksWindow(sumYspec, 100);
  % peakInd = peakInd{1};
  % peakMag = peakMag{1};
  % subplot(336);
  % plot(specInd_ica(:, 1), sumYspec, specInd_ica(:, 1)(peakInd), sumYspec(peakInd), 'or')
  % [maxV, maxI] = max(peakMag);
  % title(['sum fft ' num2str(specInd_ica(:, 1)(peakInd(maxI))) 'BPM/mag ' num2str(sumYspec(peakInd(maxI)))])

  print(picPath);
  close f;
end