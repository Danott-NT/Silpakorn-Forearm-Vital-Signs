function _run(subjectList, settingList, AC_SD, detrendFactor, sgFilter_order, sgFilter_size)
  installLibrary();
  
  clc
  close all;
  
  file = {"normal"};
  device = "LAPSUN"
  for subjects = subjectList
    for setting = settingList
      for fileNames = file
        for timeWindow = [20]
          for timeStep = [1]
            fileName = fileNames{1};
            if (setting == 3 || setting == 4)
              fileName = [fileName "_EVM_600-6-1"];
            end

            dir = "//DAXZ-SERVER/ResearchServer/ContactFreePulseMeasurement/Data/";
            length(strfind(subject, "UBFC"))
            if (length(strfind(subject, "UBFC")) == 1)
              dataSet = 2;
            else
              dataSet = 1;
            end

            if ((length(strfind(subject, "SU_led")) == 1 && subjects <= 5) )
              dir = "//DAXZ-SERVER/FileServer/ContactFreePulseMeasurement/Data/";
            end
            nPointFFTMultiple = 8;
            disp(["------" device " " subject " " fileName " " num2str(timeWindow) "-" num2str(timeStep) "------"]);
            path = [dir "/" device "/" subject];
            % verbose = 0 ==> no log
            % verbose = 1 ==> log on txt file
            % verbose = 2 ==> plot figure
            % verbose = 3 ==> log on txt and plot figure
            verbose = 3
            doICA = 1;
            if (length(strfind(fileName, "EVM")) != 0) 
              doICA = 0;
            end

            p = [dir "/" device "/" subject "/Feature/" fileName "_Skin/_feature_" fileName "_Skin.txt"]
            if (exist(p) == 0)
              continue
            end
            [flag, pulseRefArduino] = preparePulseArduino(path, fileName);
            minBPM = 45; maxBPM = 240;
            timeWindow = double(timeWindow);

            bloodVesselDamageDetection(dir, device, subject, fileName, timeWindow, timeStep, doICA, verbose, nPointFFTMultiple, minBPM, maxBPM, setting, dataSet);
            
            close all;
          end
        end
      end
    end
  end
end