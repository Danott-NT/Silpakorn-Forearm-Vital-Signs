function plotAllResultSignal(hrAllWindow, signalName, picPath, device, subject, filename, rangeBPM, timePerWindow, timeStep)
  f = figure('visible','off');
  % if (isDetrend==1)
  %   picPath = [picPath "/_bloodVessel/move" num2str(moveAVG) "-detrend/"]  
  % else
  %   picPath = [picPath "/_bloodVessel/move" num2str(moveAVG) "/"]
  % end
  if (exist(picPath) == 0)
    mkdir(picPath);
  end
  [nWindows, nSignal] = size(hrAllWindow);
  nSignal = (nSignal-2)/4; %each signal contains 4 series including magnitude, highest mag - 2nd, snr mag, and hr

  ref = hrAllWindow(:, 2); %ref magnitude is dummy (not use)
  for ind = 1:nSignal
    hr = hrAllWindow(:, 4*ind+2);
    mag = hrAllWindow(:, 4*ind+2-3);
    diffFirstSecondMag = hrAllWindow(:, 4*ind+2-2);
    snr = hrAllWindow(:, 4*ind+2-1);

    d = mean(abs(ref-hr));
    t = 1:length(ref);

    subplot(411);
    plot(ref, "marker", 'o', "markersize", 3, 'color', [1 0 0])
    hold on;
    plot(hr, "marker", 'o', "markersize", 3, 'color', [0 0 1])
    title(["HR-meanErr " num2str(d) " bpm"]);
    legend("ref", signalName{ind});

    subplot(412);
    plot(mag, "marker", 'o', "markersize", 3, "markeredgecolor", 'k', 'color', 'k'); hold on; 
    plot(diffFirstSecondMag, "marker", 'o', "markersize", 3, "markeredgecolor", 'g', 'color', 'g');
    legend('mag', '1st-2nd mag'); title('mag'); xlabel('#window')

    subplot(413); plot(snr); legend('snr'); xlabel('#window')

    subplot(427);
    d = hr-ref;
    m = (hr+ref)/2;
    plot(m, d, 'ob', 'markersize', 3);
    title('BA-hr-vs-ref')
    xlabel('mean(hr,ref)'); ylabel('hr-ref');

    subplot(428);
    plot(mag, snr, 'ok', 'markersize', 3);
    title('Scatter-mag-vs-snr')
    xlabel('mag'); ylabel('snr');
    

    % subplot(414);
    % plot(diff(hrAllWindow(:, 2*(ind+1)-1)), "marker", 'o', "markersize", 3, "markeredgecolor", 'k', 'color', [0 0 1])
    % title(["diff mag SD " num2str(std(diff(hrAllWindow(:, 2*(ind+1)-1))))]);
    % legend(signalName{ind});

    p = [picPath "/" signalName{ind} "/filter" num2str(rangeBPM) "_time" num2str(timePerWindow) "-" num2str(timeStep) "/"];
    if (exist(p) == 0)
      mkdir(p);
    end
    print([p device '_' filename "_filter" num2str(rangeBPM) "_time" num2str(timePerWindow) "-" num2str(timeStep) "_bpm_" signalName{ind} "_" subject ".png"]);
    clf;
    dlmwrite([p device '_' filename "_filter" num2str(rangeBPM) "_time" num2str(timePerWindow) "-" num2str(timeStep) "_bpm_" signalName{ind} "_" subject ".txt"]...
      , [ref, hr, mag, diffFirstSecondMag, snr], "\t");
  end
  close f;
end