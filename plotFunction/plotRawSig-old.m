function plotRawSig(windowedSig_nor, windowedSig_specIndRaw, windowedSig_ySpecRaw, peakInd, magMaxInd, pulseMaxInd, windowedSig_specIndRaw_snr, windowedSig_ySpecRaw_snr, peakInd_snr, magMaxInd_snr, pulseMaxInd_snr, picPath)
  f = figure('visible','off');
  subplot(3,3,1); plot(windowedSig_nor(:, 1), "b"); title(['Preprocess B']);
  subplot(3,3,4); plot(windowedSig_nor(:, 2), "g"); title(['Preprocess G']);
  subplot(3,3,7); plot(windowedSig_nor(:, 3), "r"); title(['Preprocess R']);

  subplot(3,3,2); plot(windowedSig_specIndRaw(:, 1), windowedSig_ySpecRaw(:, 1), "b", windowedSig_specIndRaw(:, 1)(peakInd{1}), windowedSig_ySpecRaw(:, 1)(peakInd{1}), "or")
  title([num2str(pulseMaxInd{1}) ' BPM/mag ' num2str(magMaxInd{1})])
  subplot(3,3,5); plot(windowedSig_specIndRaw(:, 1), windowedSig_ySpecRaw(:, 2), "g", windowedSig_specIndRaw(:, 1)(peakInd{2}), windowedSig_ySpecRaw(:, 2)(peakInd{2}), "or")
  title([num2str(pulseMaxInd{2}) ' BPM/mag ' num2str(magMaxInd{2})])
  subplot(3,3,8); plot(windowedSig_specIndRaw(:, 1), windowedSig_ySpecRaw(:, 3), "r", windowedSig_specIndRaw(:, 1)(peakInd{3}), windowedSig_ySpecRaw(:, 3)(peakInd{3}), "or")
  title([num2str(pulseMaxInd{3}) ' BPM/mag ' num2str(magMaxInd{3})])

  subplot(3,3,3); plot(windowedSig_specIndRaw_snr(:, 1), windowedSig_ySpecRaw_snr(:, 1), "b", windowedSig_specIndRaw_snr(:, 1)(peakInd_snr{1}), windowedSig_ySpecRaw_snr(:, 1)(peakInd_snr{1}), "or")
  title([num2str(pulseMaxInd_snr{1}) ' BPM/mag ' num2str(magMaxInd_snr{1})])
  subplot(3,3,6); plot(windowedSig_specIndRaw_snr(:, 2), windowedSig_ySpecRaw_snr(:, 2), "g", windowedSig_specIndRaw_snr(:, 1)(peakInd_snr{2}), windowedSig_ySpecRaw_snr(:, 2)(peakInd_snr{2}), "or")
  title([num2str(pulseMaxInd_snr{2}) ' BPM/mag ' num2str(magMaxInd_snr{2})])
  subplot(3,3,9); plot(windowedSig_specIndRaw_snr(:, 3), windowedSig_ySpecRaw_snr(:, 3), "r", windowedSig_specIndRaw_snr(:, 1)(peakInd_snr{3}), windowedSig_ySpecRaw_snr(:, 3)(peakInd_snr{3}), "or")
  title([num2str(pulseMaxInd_snr{3}) ' BPM/mag ' num2str(magMaxInd_snr{3})])

  print(picPath);
  close f;
end