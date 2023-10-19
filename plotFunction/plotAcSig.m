function plotAcSig(windowedSig, timePerWindow, windowedSig_peakInd, windowedSig_specIndRaw, windowedSig_ySpecRaw, peakInd, magMaxInd, pulseMaxInd, windowedSig_specIndRaw_sum, windowedSig_ySpecRaw_sum, peakInd_sum, magMaxInd_sum, pulseMaxInd_sum, picPath)
  f = figure('visible','off');
  ch = size(windowedSig, 2);

  colorPalette = {"b", "g", "r"};
  for i = [1:ch]
    
    sig = windowedSig(:, i);
    timeInd = 1:length(sig);

    subplot(3,3,3*i-2); 
    d = diff(windowedSig_peakInd{i});
    medIBI = double(median(d) * timePerWindow) / length(sig);
    pulsePeak = 60/medIBI;

    % d
    % [mean(d), mean(d)+1.5*std(d), mean(d)-1.5*std(d)]
    % plot(timeInd, sig, "b", timeInd(windowedSig_peakInd{i}), sig(windowedSig_peakInd{i}), "or"); 
    % title(['AC Com' num2str(i) ' ' num2str(length(windowedSig_peakInd{i})) ' ' num2str(pulsePeak) ' ' num2str(std(d))]);
    plot(timeInd, sig, colorPalette{i}, timeInd(windowedSig_peakInd{i}), sig(windowedSig_peakInd{i}), ".k", "markersize", 7); 
    title(['AC Com' num2str(i) ' SD:' num2str(round(std(d)*10)/10)]);
    set(gca,'XTick', floor(linspace(timeInd(1), timeInd(end), 5)));
    set(gca,'YTick', round(linspace(min(sig), max(sig), 5)*10)/10);

    subplot(3,3,3*i-1); plot(windowedSig_specIndRaw, windowedSig_ySpecRaw(:, i), "b", windowedSig_specIndRaw(peakInd{i}), windowedSig_ySpecRaw(:, i)(peakInd{i}), "or")
    title([num2str(round(pulseMaxInd{i}*10)/10) ' BPM/mag ' num2str(round(magMaxInd{i}*10)/10)], 'FontSize', 8)
    
    set(gca,'XTick', floor(linspace(windowedSig_specIndRaw(1), windowedSig_specIndRaw(end), 5)));
    set(gca,'YTick', floor(linspace(min(windowedSig_ySpecRaw(:, i)), max(windowedSig_ySpecRaw(:, i)), 5)));
  end
  % sum of fft from all comp.
  subplot(336);
  plot(windowedSig_specIndRaw_sum, windowedSig_ySpecRaw_sum, "b", windowedSig_specIndRaw_sum(peakInd_sum{1}), windowedSig_ySpecRaw_sum(peakInd_sum{1}), 'or')
  title([num2str(round(pulseMaxInd_sum{1}*10)/10) ' BPM/mag ' num2str(round(magMaxInd_sum{1}*10)/10)], 'FontSize', 8)
  set(gca,'XTick', floor(linspace(windowedSig_specIndRaw_sum(1), windowedSig_specIndRaw_sum(end), 5)));
  set(gca,'YTick', floor(linspace(min(windowedSig_ySpecRaw_sum), max(windowedSig_ySpecRaw_sum), 5)));
  % sumYspec = windowedSig_ySpecRaw(:, 1) + windowedSig_ySpecRaw(:, 2) + windowedSig_ySpecRaw(:, 3);
  % [peakMag, peakInd, pulseRefPeak, HRVpeak] = findpeaksWindow(sumYspec, 100);
  % peakInd = peakInd{1};
  % peakMag = peakMag{1};
  % subplot(336);
  % plot(windowedSig_specIndRaw(:, 1), sumYspec, windowedSig_specIndRaw(:, 1)(peakInd), sumYspec(peakInd), 'or')
  % [maxV, maxI] = max(peakMag);
  % title(['sum fft ' num2str(windowedSig_specIndRaw(:, 1)(peakInd(maxI))) 'BPM/mag ' num2str(sumYspec(peakInd(maxI)))])

  print(picPath);
  close f;
end