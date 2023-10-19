function deleterOldFolder()
	file = {"50-normal", "50-press", "180-press"...
          , "50-normal_GHX", "50-press_GHX", "180-press_GHX"};
	% ramaFile = {"normal_GHX"};
	% ramaFile = {"normal", "vein-occlusion1", "artery-occlusion"...
	%         , "normal_GHS", "vein-occlusion1_GHS", "artery-occlusion_GHS"...
	%         , "normal_GHX", "vein-occlusion1_GHX", "artery-occlusion_GHX"...
	%         , "normal_XoverY", "vein-occlusion1_XoverY", "artery-occlusion_XoverY"};
	device = "LAPSUN"
	for subjects = [1]
		for fileNames = file
		  	for timeWindow = [10]
				for timeStep = [1]
				  fileName = fileNames{1};
				  subject = ["S" num2str(subjects)];
				  nPointFFTMultiple = 3;
				  disp(["------" device " " subject " " fileName " " num2str(timeWindow) "-" num2str(timeStep) "------"]);
				  path = [dir "/" device "/" subject];
				  minBPM = 45; maxBPM = 240;
				  close all;
				end
		  	end
		end
	end
end
