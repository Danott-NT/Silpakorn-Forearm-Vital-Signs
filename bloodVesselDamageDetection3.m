function accR = bloodVesselDamageDetection(dir, device, subject, filename, timePerWindow, timeStep, doICA, verbose, nPointFFTMultiple, minBPM, maxBPM, setting, dataSet)
  
  % graphics_toolkit("gnuplot");
  graphics_toolkit("qt");
  tic % Timer
  % Path setting
  mainPath = [dir device '/' subject "/Feature/" filename "_Skin"];
  propertyPath = [mainPath '/_WindowsProperty.property'];
  [key, value] = textread(propertyPath, "%s %f");
  for i = 1:length(key)
    if strcmp(key(i), "frameRate")
      samplingRate = value(i)
    end
  end

  colorFeatures = {"BGR", "HLS", "HSV", "Lab", "Luv" , "TSL",  "YCrCb"};
  featurePath = [mainPath "/_feature_" filename "_Skin.txt"];
  % Change all dimension to ****Rows = samples, Cols = channel****
  sig = load(featurePath); %Rows = samples, Cols = channel
  channels = size(sig, 2)-1; % Minus col Frame number
  sig = sig(:, 2:channels+1); % rows : fps*60 cols : 7 colorModel * 3 Channels
  % sig = sig(1:35*60, :);
  nsamples = size(sig, 1);
  
  timeRecorded = uint16(nsamples/samplingRate);
  nsamples = uint16(timeRecorded*samplingRate);
  
  windowStep = ceil(timeStep*samplingRate);

  filterType = "Normal_Fil";
  pulsePath = [mainPath "/BloodVessel_time_" num2str(timePerWindow) "-" num2str(timeStep)];
  pulsePath = [pulsePath "_NF"];
  rangeBPM = 0;

  if (timePerWindow > timeRecorded)
    pulsePath = [pulsePath "_MAX-LENGTH"];
    timePerWindow = timeRecorded;
  end

  minHz = (minBPM-10)/60
  maxHz = (maxBPM+10)/60

  if (setting == 1) %poh+
    prefix = "poh+";
  elseif (setting == 2) 
    prefix = "poh+5";
  elseif (setting == 3)
    prefix = "wu+5";
  elseif (setting == 4)
    prefix = "wu+sg"
  elseif (setting == 5)
    prefix = "ica+move5";
  elseif (setting == 6)
    prefix = "ica+sg";
  elseif (setting == 7)
    prefix = "ica+sg_ACF";
  end
  
  pulsePath = [pulsePath "_" num2str(minBPM) "-" num2str(maxBPM) "_" prefix "/"];

  picPath = [pulsePath "pulsePic"];
  if (verbose != 0)
    if (exist(pulsePath) == 7)
      disp(["deleting old picPath..." pulsePath])
      confirm_recursive_rmdir(0)
      rmdir(picPath, 's');
      disp(["deleting old pulsePath..." pulsePath])
      % [status, msg, msgid] = rmdir(pulsePath, 's')
      pause(2);
    end
  end
  disp(["mkdir pulsePath..." pulsePath]);
  [STATUS, MSG] = mkdir(pulsePath)
  
  disp(["mkdir picPath..." pulsePath]);
  mkdir(picPath)
  

  arduinoSync = load([dir device '/' subject "/Video_Frame/" filename "_ori/_rawPulseArduino_Sync.txt"]);
%  Arduino pulse format : column 1:hr, 2:min, 3:sec, 4:signal(use), 5:signal, 6:hr
  arduinoSync = fixArduinoFPS(arduinoSync);
  timeSec = arduinoSync(:, 1);
  
  timePerWindow = min([timePerWindow, timeSec(end), timeRecorded]);
  samplePerWindow = timePerWindow*samplingRate;
  windowedSig = zeros(samplePerWindow, 3);
  nWindowsVideo = uint16((timeRecorded-timePerWindow+1)/timeStep)
  nWindowsSignal = (size(sig, 1)-timePerWindow*samplingRate+1)/(timeStep*samplingRate)
  nWindowsArduino = uint16((timeSec(end)-timePerWindow+1)/timeStep)-1
  nWindows = min([nWindowsVideo, nWindowsArduino, nWindowsSignal])
  totalTime = nWindows * timeStep;
  
  hrAllWindow = [];

  % nWindows = 20;
  currPlot = 1;
  nPlot = nWindows;
  if (nPlot > 40)
    nPlot = 40;
  end
  
  hrAllWindow = [];
  colorModel = "BGR";
  % Sliding window without loop can be acheived by using im2col or nlfilter
  % source : https://stackoverflow.com/questions/5262630/matlab-moving-window-avoiding-nested-loops
  sig = sig'; %row = ch, col = sample
  sig = sig(1:3, :); %select only RGB
  dt = 1/samplingRate;
  t = (0:length(sig)-1)/samplingRate;
  % samplingRate = 256;
  dt2 = 1/samplingRate;
  % samplePerWindow = timePerWindow*samplingRate;
  t2 = 0:dt2:t(end)-dt2;
  % sig = [interp1(t, sig(1, :), t2, "spline"); interp1(t, sig(2, :), t2, "spline"); interp1(t, sig(3, :), t2, "spline")];

  sig_time = double(1:size(sig, 2))/double(samplingRate);
  % dlmwrite("C:/Users/DAXZ-GAMING/Desktop/research/reportAC/1.txt", sig', "\t")
  % Sliding window on whole data
  % im2col explanation : http://solvecode.blogspot.com/2016/06/img2col-function-in-matlab-explanation.html
  pkg load image
  % x = 1:length(sig);
  % samplingRate = 64;
  % t2 = 1/samplingRate;
  % x2 = 1:t2:sig_time(end);
  % sig = cell2mat(arrayfun(@(i) interp1(sig_time, sig(i, :), x2), 1:size(sig, 1), 'UniformOutput', false)');
  % sig_time = x2;

  % ******************************************************
  % Windowing RAW RGB ** USE FOR PLOT**
  windowedSig_raw = im2col(sig, [1, samplePerWindow], 'sliding'); %row = sample, col = ch
  % Split slided data into cell, each cell contain 3channel of sliding window
  windowedSig_raw = mat2cell(windowedSig_raw, [samplePerWindow], 3*ones(1, size(windowedSig_raw, 2)/3)); %row = channel, col = ch
  windowedSig_raw = windowedSig_raw(1:windowStep:end);
  if (nWindows > length(windowedSig_raw))
    nWindows = length(windowedSig_raw);
  end
  windowedSig_raw = windowedSig_raw(1:nWindows);
  % ******************************************************

  windowedSig = im2col(sig, [1, samplePerWindow], 'sliding'); %row = sample, col = ch
  % Split slided data into cell, each cell contain 3channel of sliding window
  windowedSig = mat2cell(windowedSig, [samplePerWindow], 3*ones(1, size(windowedSig, 2)/3)); %row = channel, col = ch
  windowedSig = windowedSig(1:windowStep:end);
  if (nWindows > length(windowedSig))
    nWindows = length(windowedSig);
  end
  windowedSig = windowedSig(1:nWindows);
  % Create time index
  windowedSig_time = im2col(sig_time, [1, samplePerWindow], 'sliding');
  windowedSig_time = mat2cell(windowedSig_time, [samplePerWindow], ones(1, size(windowedSig_time, 2)));
  windowedSig_time = windowedSig_time(1:windowStep:end);
  windowedSig_time = windowedSig_time(1:nWindows);
  
  disp(["timePerWindow : " num2str(timePerWindow) " nWindows : " num2str(nWindows) " nPlot :" num2str(nPlot)])
  %%%%%%%%%%%%%%%%%%%           Arduino Ref
  % timeStart = 1;
  [arStart, arStop] = arrayfun(@(timeStart) splitArduinoIndexToWindow(timeSec, timeStart, timeStep, timePerWindow), 1:timeStep:totalTime);
  arStart = arStart(1:nWindows);
  arStop = arStop(1:nWindows);
  pulse_fps_window = ((arStop-arStart)+1)/double(timePerWindow);
  a = prepareData(arduinoSync(:, 4), mean(pulse_fps_window), 'sgFilter', 2, 127, 'filter', [minBPM, maxBPM]);
  arduinoSync(:, 4) = a;
  arduinoSignal = arrayfun(@(startInd, stopInd, p_fps) prepareData(arduinoSync(startInd:stopInd, 4), p_fps, 'normalize'), arStart, arStop, pulse_fps_window, 'UniformOutput', false);
  pulseRefArduinoRaw = arrayfun(@(startInd, stopInd) arduinoSync(startInd:stopInd, 6), arStart, arStop, 'UniformOutput', false);
  pulseRef_arduino = cellfun(@(p) mean(p), pulseRefArduinoRaw);

  pulseRefSD_arduino = arrayfun(@(startInd, stopInd) std(arduinoSync(startInd:stopInd, 6)), arStart, arStop);
  [pulseRefFFT_arduino, pulseRefFFT_mag, pulseRefFFT_bpm] = cellfun(@(x, p_fps) specAnalysis(x, p_fps, minBPM, maxBPM, nPointFFTMultiple), arduinoSignal, num2cell(pulse_fps_window), 'UniformOutput', false);
  [peakMag_arduino, peakInd_arduino, pulseRefPeak_arduino, HRVpeak_arduino] = cellfun(@(x, p_fps) findpeaksWindow(x, p_fps, 'timeStep', 0.2), arduinoSignal, num2cell(pulse_fps_window), 'UniformOutput', false);
  % [peakMag_arduino, peakInd_arduino, pulseRefPeak_arduino, HRVpeak_arduino] = cellfun(@(x) findpeaksWindow(x, timePerWindow), arduinoSignal, 'UniformOutput', false);
 
  checkWindowInd = [];
  pulseRef = [];
  pulseRef = pulseRefPeak_arduino{1}{1};
  % pulseRef = pulseRefFFT_arduino{1}(1, 2)
  lastRef = pulseRef;
  % pulseCheckCount = 0;
  for ind = 2:nWindows
    p_ar = pulseRef_arduino(ind);
    p_fft = pulseRefFFT_arduino{ind}(1, 2);
    p_peak = pulseRefPeak_arduino{ind}{1};

    if (length(pulseRef) > 5)
      lastRef = median(pulseRef(end-5:end));
    end
    
    if (abs(p_fft - p_peak) < 5)
      pulseRef_tmp = p_fft;
    elseif (abs(lastRef - p_peak) < 5)
      pulseRef_tmp = p_peak;
    elseif (abs(lastRef - p_fft) < 5)
      pulseRef_tmp = p_fft;
    else
      pulseRef_tmp = lastRef;
      checkWindowInd = [checkWindowInd, ind];
      if (ind-1 >= 1)
        checkWindowInd = [checkWindowInd, ind-1];
      end
      if (ind+1 <= nWindows)
        checkWindowInd = [checkWindowInd, ind+1];
      end
    end
    pulseRef = [pulseRef, pulseRef_tmp];
  end
  checkWindowInd = unique(checkWindowInd);
  pulseRef

  if (verbose >= 2)
    arrayfun(@(ind) plotArduinoSig_fixedAxes(ind, arduinoSignal{ind}, pulse_fps_window(ind), timePerWindow...
        , pulseRefArduinoRaw{ind}, pulseRef_arduino(ind) ...
        , pulseRefFFT_arduino{ind}(1, 1:2), pulseRefFFT_mag{ind}, pulseRefFFT_bpm{ind}...
        , pulseRefPeak_arduino{ind}, peakInd_arduino{ind}...
        , [picPath "/" num2str(ind) "_signal-arduino.png"]), [1:nPlot, checkWindowInd]);
  end
  checkWindowInd(checkWindowInd==0) = [];

  % Set filter range for "CUSTOM" parameter
  if strcmp (filterType, "Custom_Fil")
    minBPM = arrayfun(@(x) x-rangeBPM, pulseRef);
    maxBPM = arrayfun(@(x) x+rangeBPM, pulseRef);
  else
    minBPM = ones(1, nWindows)*minBPM;
    maxBPM = ones(1, nWindows)*maxBPM;
  end
  minBPM = mat2cell(minBPM, 1, ones(1, nWindows));
  maxBPM = mat2cell(maxBPM, 1, ones(1, nWindows));

  %%%%%%%%%%%%%%%%%%%%%%%

  % Do Function on each column of cell array
  % cellfun(function, data) ==> apply function on each cell
  if (setting == 1) %poh+
    disp("poh+")
    windowedSig_nor = cellfun(@(x, minHR, maxHR) prepareData(x, samplingRate, 'detrendLambda', 10, 'normalize'), windowedSig, minBPM, maxBPM, 'UniformOutput',false);
  elseif (setting == 2) %ica
    windowedSig_nor = cellfun(@(x, minHR, maxHR) prepareData(x, samplingRate, 'detrendLambda', 20, 'normalize'), windowedSig, minBPM, maxBPM, 'UniformOutput',false);
  elseif (setting == 3) %EVM
    disp("wu+")
    windowedSig_nor = cellfun(@(x, minHR, maxHR) prepareData(x, samplingRate, 'moveAVG', 5, 'filter', [minHR, maxHR], 'normalize'), windowedSig, minBPM, maxBPM, 'UniformOutput',false);
  elseif (setting == 4) %EVM sg
    disp("wu+sg")
    windowedSig_nor = cellfun(@(x, minHR, maxHR) prepareData(x, samplingRate, 'sgFilter', 4, 25, 'filter', [minHR, maxHR], 'normalize'), windowedSig, minBPM, maxBPM, 'UniformOutput',false);
  elseif (setting >= 5) % OUR Method
    disp("ica+")
    if (dataSet == 1) %our
      lambda = 20
    else  %UBFC
      lambda = 35
    end
    windowedSig_nor = cellfun(@(x, minHR, maxHR) prepareData(x, samplingRate, 'detrendLambda', lambda, 'normalize'), windowedSig, minBPM, maxBPM, 'UniformOutput',false);

  % elseif (setting == 4) %ica+
  %   if (dataSet == 1) %our
  %     disp("ica+ dataSet1")
  %     windowedSig_nor = cellfun(@(x, minHR, maxHR) prepareData(x, samplingRate, 'detrendLambda', 20, 'normalize'), windowedSig, minBPM, maxBPM, 'UniformOutput',false);
  %   elseif (dataSet == 2) %UBFC
  %     disp("ica+ dataSet2")
  %     windowedSig_nor = cellfun(@(x, minHR, maxHR) prepareData(x, samplingRate, 'detrendLambda', 35, 'normalize'), windowedSig, minBPM, maxBPM, 'UniformOutput',false);
  %   end
  end

  if (verbose == 1 || verbose == 3)
    % cellfun(@(x, i) dlmwrite([pulsePath num2str(i) "_" colorModel "_windowed_Signal_tab.txt"], x, "\t"), windowedSig, num2cell(1:nWindows,1));
    cellfun(@(x, i) writeTxtFile([pulsePath num2str(i) "_" colorModel "_windowed_Signal_tab.txt"], x), windowedSig, num2cell(1:nWindows,1));
    % Write data before normalization
    cellfun(@(x, i) writeTxtFile([pulsePath num2str(i) "_" colorModel "_windowed_Signal_PreProcess_tab.txt"], x), windowedSig_nor, num2cell(1:nWindows,1));
  end
  % Result is in form of nested cell. Each cell contain 3 subcell which is 3 output of specAnalysis
  [pulse_windowedSig, ySpec_windowedSig, specInd_windowedSig] = cellfun(...
        @(x, minHR, maxHR) specAnalysis(x, samplingRate, minHR, maxHR, nPointFFTMultiple)...
        , windowedSig_nor, minBPM, maxBPM, 'UniformOutput',false);
  [peakMag_windowedSig, peakInd_windowedSig, pulsePeak_windowedSig, HRVpeak_windowedSig] = cellfun(...
        @(x) findpeaksWindow(x, samplingRate, 'timeStep', 0.2)...
        , ySpec_windowedSig, 'UniformOutput',false);
  % Print pulse
  if (verbose == 1 || verbose == 3)
    cellfun(@(x, i) writeTxtFile([pulsePath num2str(i) "_" colorModel "_windowed_Signal_pulse_tab.txt"], x)...
          , pulse_windowedSig, num2cell(1:nWindows,1));
  end
  % Find max magnitude
  peakMaxInd_windowedSig = cellfun(@(magCell, indCell) ...
        cellfun(@(mag, ind) ind(find(mag == max(mag))), magCell, indCell, 'UniformOutput',false)...
        , peakMag_windowedSig, peakInd_windowedSig, 'UniformOutput',false);
  [~, sortI_windowedSig] = cellfun(@(magCell) ...
        cellfun(@(mag) sort(mag, 'descend'), magCell, 'UniformOutput',false)...
        , peakMag_windowedSig, 'UniformOutput',false);
  firstPeakMag_windowedSig = cellfun(@(magCell, indCell, sortICell) ...
        cellfun(@(mag, ind, sortI) mag(sortI(1)), magCell, indCell, sortICell, 'UniformOutput',false)...
        , peakMag_windowedSig, peakInd_windowedSig, sortI_windowedSig, 'UniformOutput',false);
  
  secondPeakMag_windowedSig = cellfun(@(magCell, indCell, sortICell) ...
        cellfun(@(mag, ind, sortI) mag(sortI(2)), magCell, indCell, sortICell, 'UniformOutput',false)...
        , peakMag_windowedSig, peakInd_windowedSig, sortI_windowedSig, 'UniformOutput',false);
  diffFirstSecodnPeakMag_windowedSig = cellfun(@(firstPeakMagCell, secondPeakMagCell) ...
        cellfun(@(f, s) f-s, firstPeakMagCell, secondPeakMagCell, 'UniformOutput',false)...
        , firstPeakMag_windowedSig, secondPeakMag_windowedSig, 'UniformOutput',false);
  snr_windowedSig = cellfun(@(mag_array, specInd, peakInd_cell) ...
        calSNRMaxPeak(mag_array, specInd, cell2mat(peakInd_cell))...
        , ySpec_windowedSig, specInd_windowedSig, peakMaxInd_windowedSig, 'UniformOutput',false);
  
  [magMaxInd_windowedSig, pulseMaxInd_windowedSig] = cellfun(@(mag_array, specInd, peakInd_cell)... 
    arrayfun(@(ch) deal(mag_array(peakInd_cell{ch}, ch), specInd(peakInd_cell{ch})), 1:size(peakInd_cell, 2), 'UniformOutput',false)...
    , ySpec_windowedSig, specInd_windowedSig, peakMaxInd_windowedSig, 'UniformOutput',false);
  
  % Sum of mag ============== test
  ySpec_sum_windowedSig = cellfun(@(x) sum(x, 2), ySpec_windowedSig, 'UniformOutput',false);
  % Print FFT
  if (verbose == 1 || verbose == 3)
    cellfun(@(ySpec, ySpecSum, specInd, i) writeTxtFile([pulsePath num2str(i) "_" colorModel "_windowed_Signal_FFT_tab.txt"], [specInd, ySpec, ySpecSum])...
          , ySpec_windowedSig, ySpec_sum_windowedSig, specInd_windowedSig, num2cell(1:nWindows,1));
  end
  specInd_sum_windowedSig = specInd_windowedSig;
  [peakMag_sum_windowedSig, peakInd_sum_windowedSig, ~, ~] = cellfun(@(x) findpeaksWindow(x, samplingRate, 'timeStep', 0.2)...
        , ySpec_sum_windowedSig, 'UniformOutput',false);
  peakMaxInd_sum_windowedSig = cellfun(@(magCell, indCell) ...
        cellfun(@(mag, ind) ind(find(mag == max(mag))), magCell, indCell, 'UniformOutput',false)...
        , peakMag_sum_windowedSig, peakInd_sum_windowedSig, 'UniformOutput',false);
  [~, sortI_sum_windowedSig] = cellfun(@(magCell) ...
        cellfun(@(mag) sort(mag, 'descend'), magCell, 'UniformOutput',false)...
        , peakMag_sum_windowedSig, 'UniformOutput',false);
  firstPeakMag_sum_windowedSig = cellfun(@(magCell, indCell, sortICell) ...
        cellfun(@(mag, ind, sortI) mag(sortI(1)), magCell, indCell, sortICell, 'UniformOutput',false)...
        , peakMag_sum_windowedSig, peakInd_sum_windowedSig, sortI_sum_windowedSig, 'UniformOutput',false);
  secondPeakMag_sum_windowedSig = cellfun(@(magCell, indCell, sortICell) ...
        cellfun(@(mag, ind, sortI) mag(sortI(2)), magCell, indCell, sortICell, 'UniformOutput',false)...
        , peakMag_sum_windowedSig, peakInd_sum_windowedSig, sortI_sum_windowedSig, 'UniformOutput',false);
  diffFirstSecodnPeakMag_sum_windowedSig = cellfun(@(firstPeakMagCell, secondPeakMagCell) ...
        cellfun(@(f, s) f-s, firstPeakMagCell, secondPeakMagCell, 'UniformOutput',false)...
        , firstPeakMag_sum_windowedSig, secondPeakMag_sum_windowedSig, 'UniformOutput',false);
  snr_sum_windowedSig = cellfun(@(mag_array, specInd, peakInd_cell) ...
        calSNRMaxPeak(mag_array, specInd, cell2mat(peakInd_cell))...
        , ySpec_sum_windowedSig, specInd_sum_windowedSig, peakMaxInd_sum_windowedSig, 'UniformOutput',false);
  [magMaxInd_sum_windowedSig, pulseMaxInd_sum_windowedSig] = cellfun(@(mag_array, specInd, peakInd_cell)... 
    arrayfun(@(ch) deal(mag_array(peakInd_cell{ch}, ch), specInd(peakInd_cell{ch})), 1:size(peakInd_cell, 2), 'UniformOutput',false)...
    , ySpec_sum_windowedSig, specInd_sum_windowedSig, peakMaxInd_sum_windowedSig, 'UniformOutput',false);
  

  if (verbose >= 2)
    [peakMag_raw, peakInd_raw, pulseRefPeak_raw, HRVpeak_raw] = cellfun(@(x) findpeaksWindow(x, samplingRate, 'timeStep', 0.2 ), windowedSig_nor, 'UniformOutput', false);
    arrayfun(@(ind) plotArduinoSig(ind, pulseRef(ind), arduinoSignal{ind}...
      , pulseRefArduinoRaw{ind}, pulseRef_arduino(ind), pulseRefSD_arduino(ind), pulseRefFFT_arduino{ind}(1, 1:2)...
      , pulseRefPeak_arduino{ind}, peakInd_arduino{ind}, windowedSig_nor{ind}, windowedSig_time{ind}, peakInd_raw{ind}...
      , specInd_windowedSig{ind}, ySpec_windowedSig{ind}, peakInd_windowedSig{ind}, peakMaxInd_windowedSig{ind}...
      , specInd_windowedSig{ind}, ySpec_sum_windowedSig{ind}, peakInd_sum_windowedSig{ind}, peakMaxInd_sum_windowedSig{ind}...
      , [picPath "/" num2str(ind) "_signal-rgb.png"]), [1:nPlot, checkWindowInd]);
    arrayfun(@(ind) plotRawSig(ind, pulseRef(ind), arduinoSignal{ind}, peakInd_arduino{ind}, windowedSig_raw{ind}, windowedSig_nor{ind}, windowedSig_time{ind}...
      , [picPath "/" num2str(ind) "_signal-raw.png"]), [1:nPlot, checkWindowInd]);
  end
  % AutoCorr
  disp("raw AC")
  % Creat ac signal
  [autocor, lags] = cellfun(@(x) ...
        arrayfun(@(ind) xcorr(x(:, ind)), 1:size(x, 2), 'UniformOutput', false)...
        , windowedSig_nor, 'UniformOutput', false);
  % Export ac signal
  % autocor_print = cellfun(@(cor_cell) cellfun(@(cor) ((cor-min(cor))/(max(cor)-min(cor))*2-1), cor_cell, 'UniformOutput', false), autocor, 'UniformOutput', false);
  autocor = cellfun(@(cor_cell) cell2mat(cor_cell), autocor, 'UniformOutput', false);
  lags = cellfun(@(l_cell) cellfun(@(l) l', l_cell, 'UniformOutput', false), lags, 'UniformOutput', false);
  lags = cellfun(@(l_cell) cell2mat(l_cell), lags, 'UniformOutput', false);
  if (verbose == 1 || verbose == 3)
    cellfun(@(x, i) writeTxtFile([pulsePath num2str(i) "_" colorModel "_windowed_Signal_AC.txt"], x), autocor, num2cell(1:nWindows,1));
  end
  % Select half of AC
  autocor = cellfun(@(x) x(int16(end/2):end, :), autocor, 'UniformOutput', false);
  lags = cellfun(@(x) x(int16(end/2):end, :), lags, 'UniformOutput', false);
  
  [pulse_windowedSig_ac, ySpec_windowedSig_ac, specInd_windowedSig_ac] = cellfun(@(cor, minHR, maxHR) specAnalysis(cor, samplingRate, minHR, maxHR, nPointFFTMultiple), autocor, minBPM, maxBPM, 'UniformOutput',false);
  [peakMag_ac, peakInd_ac, pulseRefPeak_ac, HRVpeak_ac] = cellfun(@(x) findpeaksWindow(x, samplingRate, 'timeStep', 0.2 ), autocor, 'UniformOutput', false);
  % AutoCorr Peak FFT
  [peakMag_windowedSig_ac, peakInd_windowedSig_ac, pulsePeak_windowedSig_ac, HRVpeak_windowedSig_ac] = ...
        cellfun(@(x) findpeaksWindow(x, samplingRate, 'timeStep', 0.2), ySpec_windowedSig_ac, 'UniformOutput',false);
  peakMaxInd_windowedSig_ac = cellfun(@(magCell, indCell) ...
        cellfun(@(mag, ind) ind(find(mag == max(mag))), magCell, indCell, 'UniformOutput',false)...
        , peakMag_windowedSig_ac, peakInd_windowedSig_ac, 'UniformOutput',false);
  [~, sortI_windowedSig_ac] = cellfun(@(magCell) ...
        cellfun(@(mag) sort(mag, 'descend'), magCell, 'UniformOutput',false)...
        , peakMag_windowedSig_ac, 'UniformOutput',false);
  firstPeakMag_windowedSig_ac = cellfun(@(magCell, indCell, sortICell) ...
        cellfun(@(mag, ind, sortI) mag(sortI(1)), magCell, indCell, sortICell, 'UniformOutput',false)...
        , peakMag_windowedSig_ac, peakInd_windowedSig_ac, sortI_windowedSig_ac, 'UniformOutput',false);
  secondPeakMag_windowedSig_ac = cellfun(@(magCell, indCell, sortICell) ...
        cellfun(@(mag, ind, sortI) mag(sortI(2)), magCell, indCell, sortICell, 'UniformOutput',false)...
        , peakMag_windowedSig_ac, peakInd_windowedSig_ac, sortI_windowedSig_ac, 'UniformOutput',false);
  diffFirstSecodnPeakMag_windowedSig_ac = cellfun(@(firstPeakMagCell, secondPeakMagCell) ...
        cellfun(@(f, s) f-s, firstPeakMagCell, secondPeakMagCell, 'UniformOutput',false)...
        , firstPeakMag_windowedSig_ac, secondPeakMag_windowedSig_ac, 'UniformOutput',false);
  snr_windowedSig_ac = cellfun(@(mag_array, specInd, peakInd_cell) ...
        calSNRMaxPeak(mag_array, specInd, cell2mat(peakInd_cell))...
        , ySpec_windowedSig_ac, specInd_windowedSig_ac, peakMaxInd_windowedSig_ac, 'UniformOutput',false);
  [magMaxInd_windowedSig_ac, pulseMaxInd_windowedSig_ac] = cellfun(@(mag_array, specInd, peakInd_cell) ...
        arrayfun(@(ch) deal(mag_array(peakInd_cell{ch}, ch), specInd(peakInd_cell{ch})), 1:size(peakInd_cell, 2), 'UniformOutput',false)...
        , ySpec_windowedSig_ac, specInd_windowedSig_ac, peakMaxInd_windowedSig_ac, 'UniformOutput',false);

  % Sum of AC mag
  ySpec_sum_windowedSig_ac = cellfun(@(x) sum(x, 2), ySpec_windowedSig_ac, 'UniformOutput',false);
  specInd_sum_windowedSig_ac = specInd_windowedSig_ac;
  [peakMag_sum_windowedSig_ac, peakInd_sum_windowedSig_ac, ~, ~] = cellfun(@(x) findpeaksWindow(x, samplingRate, 'timeStep', 0.2)...
        , ySpec_sum_windowedSig_ac, 'UniformOutput',false);
  peakMaxInd_sum_windowedSig_ac = cellfun(@(magCell, indCell) ...
        cellfun(@(mag, ind) ind(find(mag == max(mag))), magCell, indCell, 'UniformOutput',false)...
        , peakMag_sum_windowedSig_ac, peakInd_sum_windowedSig_ac, 'UniformOutput',false);
  [~, sortI_sum_windowedSig_ac] = cellfun(@(magCell) ...
        cellfun(@(mag) sort(mag, 'descend'), magCell, 'UniformOutput',false)...
        , peakMag_sum_windowedSig_ac, 'UniformOutput',false);
  firstPeakMag_sum_windowedSig_ac = cellfun(@(magCell, indCell, sortICell) ...
        cellfun(@(mag, ind, sortI) mag(sortI(1)), magCell, indCell, sortICell, 'UniformOutput',false)...
        , peakMag_sum_windowedSig_ac, peakInd_sum_windowedSig_ac, sortI_sum_windowedSig_ac, 'UniformOutput',false);
  secondPeakMag_sum_windowedSig_ac = cellfun(@(magCell, indCell, sortICell) ...
        cellfun(@(mag, ind, sortI) mag(sortI(2)), magCell, indCell, sortICell, 'UniformOutput',false)...
        , peakMag_sum_windowedSig_ac, peakInd_sum_windowedSig_ac, sortI_sum_windowedSig_ac, 'UniformOutput',false);
  diffFirstSecodnPeakMag_sum_windowedSig_ac = cellfun(@(firstPeakMagCell, secondPeakMagCell) ...
        cellfun(@(f, s) f-s, firstPeakMagCell, secondPeakMagCell, 'UniformOutput',false)...
        , firstPeakMag_sum_windowedSig_ac, secondPeakMag_sum_windowedSig_ac, 'UniformOutput',false);
  snr_sum_windowedSig_ac = cellfun(@(mag_array, specInd, peakInd_cell) ...
        calSNRMaxPeak(mag_array, specInd, cell2mat(peakInd_cell))...
        , ySpec_sum_windowedSig_ac, specInd_sum_windowedSig_ac, peakMaxInd_sum_windowedSig_ac, 'UniformOutput',false);
  [magMaxInd_sum_windowedSig_ac, pulseMaxInd_sum_windowedSig_ac] = cellfun(@(mag_array, specInd, peakInd_cell) ...
        arrayfun(@(ch) deal(mag_array(peakInd_cell{ch}, ch), specInd(peakInd_cell{ch})), 1:size(peakInd_cell, 2), 'UniformOutput',false)...
        , ySpec_sum_windowedSig_ac, specInd_sum_windowedSig_ac, peakMaxInd_sum_windowedSig_ac, 'UniformOutput',false);
  
  % Mag history raw
  magHistory = [pulseRef', cell2mat(cellfun(@(x) x(1), pulseMaxInd_windowedSig))', cell2mat(cellfun(@(x) x(1), magMaxInd_windowedSig))'];
  p = [dir "_magHistory/filter" num2str(rangeBPM) "_time" num2str(timePerWindow) "-" num2str(timeStep) "/"];
  if (exist(p) == 0)
    mkdirs(p);
  end
  writeTxtFile([p device '_' filename "_" subject ".txt"], magHistory);
  f = figure('visible','off');
  subplot(211); hold on;
  plot(magHistory(:, 1), 'r'); 
  plot(magHistory(:, 2), 'b'); title('HR');
  legend('ref', 'pred');

  subplot(212)
  plot(magHistory(:, 3)); title("FFT magnitude");
  print([p device '_' filename "_" subject ".png"]);
  clf;
  % Mag history ac
  magHistory = [pulseRef', cell2mat(cellfun(@(x) x(1), pulseMaxInd_windowedSig_ac))', cell2mat(cellfun(@(x) x(1), magMaxInd_windowedSig_ac))'];
  writeTxtFile([p device '_' filename "_" subject "_ac.txt"], magHistory);
  f = figure('visible','off');
  subplot(211); hold on;
  plot(magHistory(:, 1), 'r'); 
  plot(magHistory(:, 2), 'b'); title('HR');
  legend('ref', 'pred');

  subplot(212)
  plot(magHistory(:, 3)); title("FFT magnitude");
  print([p device '_' filename "_" subject "_ac.png"]);
  clf;

  if (verbose >= 2)
    arrayfun(@(ind) plotAcSig(autocor{ind}, timePerWindow, peakInd_ac{ind}...
      , specInd_windowedSig_ac{ind}, ySpec_windowedSig_ac{ind}, peakMaxInd_windowedSig_ac{ind}, magMaxInd_windowedSig_ac{ind}, pulseMaxInd_windowedSig_ac{ind} ...
      , specInd_sum_windowedSig_ac{ind}, ySpec_sum_windowedSig_ac{ind}, peakMaxInd_sum_windowedSig_ac{ind}, magMaxInd_sum_windowedSig_ac{ind}, pulseMaxInd_sum_windowedSig_ac{ind} ...
      , [picPath "/" num2str(ind) "_signal-raw_ac.png"]), [1:nPlot, checkWindowInd]);
  end

  % Create output
  % Max mag among all component windowedSig
  [maxV, maxI] = cellfun(@(mag) max([mag{:}]), magMaxInd_windowedSig, 'UniformOutput',false);
  mag_maxMag_windowedSig = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, magMaxInd_windowedSig, 'UniformOutput',false);
  pulse_maxMag_windowedSig = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, pulseMaxInd_windowedSig, 'UniformOutput',false);
  diffFirstSecond_maxMag_windowedSig = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, diffFirstSecodnPeakMag_windowedSig, 'UniformOutput',false);
  snr_maxMag_windowedSig = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, snr_windowedSig, 'UniformOutput',false);

  % Sum magnitude
  [maxV, maxI] = cellfun(@(mag) max([mag{:}]), magMaxInd_sum_windowedSig, 'UniformOutput',false);
  mag_maxMag_sum_windowedSig = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, magMaxInd_sum_windowedSig, 'UniformOutput',false);
  pulse_maxMag_sum_windowedSig = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, pulseMaxInd_sum_windowedSig, 'UniformOutput',false);
  diffFirstSecond_maxMag_sum_windowedSig = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, diffFirstSecodnPeakMag_sum_windowedSig, 'UniformOutput',false);
  snr_maxMag_sum_windowedSig = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, snr_sum_windowedSig, 'UniformOutput',false);

  % Max mag among all component windowedSig ac
  [maxV, maxI] = cellfun(@(mag) max([mag{:}]), magMaxInd_windowedSig_ac, 'UniformOutput',false);
  mag_maxMag_windowedSig_ac = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, magMaxInd_windowedSig_ac, 'UniformOutput',false);
  pulse_maxMag_windowedSig_ac = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, pulseMaxInd_windowedSig_ac, 'UniformOutput',false);
  diffFirstSecond_maxMag_windowedSig_ac = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, diffFirstSecodnPeakMag_windowedSig_ac, 'UniformOutput',false);
  snr_maxMag_windowedSig_ac = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, snr_windowedSig_ac, 'UniformOutput',false);

  % Sum AC magnitude
  [maxV, maxI] = cellfun(@(mag) max([mag{:}]), magMaxInd_sum_windowedSig_ac, 'UniformOutput',false);
  mag_maxMag_sum_windowedSig_ac = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, magMaxInd_sum_windowedSig_ac, 'UniformOutput',false);
  pulse_maxMag_sum_windowedSig_ac = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, pulseMaxInd_sum_windowedSig_ac, 'UniformOutput',false);
  diffFirstSecond_maxMag_sum_windowedSig_ac = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, diffFirstSecodnPeakMag_sum_windowedSig_ac, 'UniformOutput',false);
  snr_maxMag_sum_windowedSig_ac = cellfun(@(maxInd, x) [x{:}](maxInd), maxI, snr_sum_windowedSig_ac, 'UniformOutput',false);

  hrAllWindow = arrayfun(@(ind) ...
    [999, pulseRef(ind), mag_maxMag_windowedSig{ind}, diffFirstSecond_maxMag_windowedSig{ind}, snr_maxMag_windowedSig{ind}, pulse_maxMag_windowedSig{ind}...
    , mag_maxMag_sum_windowedSig{ind}, diffFirstSecond_maxMag_sum_windowedSig{ind}, snr_maxMag_sum_windowedSig{ind}, pulse_maxMag_sum_windowedSig{ind}...
    , mag_maxMag_windowedSig_ac{ind}, diffFirstSecond_maxMag_windowedSig_ac{ind}, snr_maxMag_windowedSig_ac{ind}, pulse_maxMag_windowedSig_ac{ind} ...
    , mag_maxMag_sum_windowedSig_ac{ind}, diffFirstSecond_maxMag_sum_windowedSig_ac{ind}, snr_maxMag_sum_windowedSig_ac{ind}, pulse_maxMag_sum_windowedSig_ac{ind}] ...
    , 1:nWindows, 'UniformOutput',false);
  
  plotName = {'maxMagRaw', 'maxMagRaw-sum', "maxMagRaw-ac", "maxMagRaw-ac-sum"};

  %%%%%%%%%%%%%%%% ICA signal - jadeR
  if (doICA)
    [hr, report_ac] = createSigICA2(arduinoSignal, windowedSig_nor, pulseRef, minBPM, maxBPM, samplingRate, timePerWindow, nPointFFTMultiple, pulsePath, colorModel, picPath, "jadeR", checkWindowInd, setting, verbose=verbose, nPlot=nPlot);
        
    hrAllWindow = cellfun(@(x, y) [x, y] , hrAllWindow, hr, 'UniformOutput',false);
    plotName = [plotName, {'maxMag-jadeR', 'maxMag-jadeR-sum', 'maxMag-jadeR-ac', 'maxMag-jadeR-ac-sum'}];
  end
  [~, nResult] = size(hrAllWindow{1});
  hrAllWindow = cell2mat(hrAllWindow);
  hrAllWindow = reshape(hrAllWindow, nResult, size(hrAllWindow, 2)/nResult)';
  writeTxtFile([pulsePath "__HRAllwindow.txt"], hrAllWindow(:, 2:4:end));
  writeTxtFile([pulsePath "__HRAllwindow_with_magnitude.txt"], hrAllWindow);

  
  if (setting == 1) %poh+
    prefix = "poh+";
  elseif (setting == 2) 
    prefix = "poh+5";
  elseif (setting == 3)
    prefix = "wu+5";
  elseif (setting == 4)
    prefix = "wu+sg";
  elseif (setting == 5)
    prefix = "ica+move5";
  elseif (setting == 6)
    prefix = "ica+sg";
  elseif (setting == 7)
    prefix = "ica+sg_ACF";
  end

  picPath = [dir "/_bloodVessel/" prefix]
  plotAllResultSignal(hrAllWindow, plotName, picPath, device, subject, filename, rangeBPM, timePerWindow, timeStep);
  
  f = figure('visible','off');
  ref = hrAllWindow(:, 2);
  nPlot = (size(hrAllWindow, 2)-2)/4;
  for ind = 1:nPlot
    hr = hrAllWindow(:, 4*ind+2);
    mag = hrAllWindow(:, 4*ind+2-3);
    diffFirstSecondMag = hrAllWindow(:, 4*ind+2-2);
    snr = hrAllWindow(:, 4*ind+2-1);
    subplot(311); plot(ref, 'r'); hold on; plot(hr); 
    title(['hr-meanErr ' num2str(mean(abs(hr-ref))) ' bpm']); legend('ref', plotName{ind}); 
    subplot(312); plot(mag); hold on; plot(diffFirstSecondMag); legend('mag', '1st-2nd mag'); title('mag'); xlabel('#window')
    subplot(313); plot(snr); legend('snr'); xlabel('#window')
    print([pulsePath "__" plotName{ind} ".png"]);
    clf;
  end
  
  totalTime_sec = toc
  totalTime_min = totalTime_sec/60
  close f;
end