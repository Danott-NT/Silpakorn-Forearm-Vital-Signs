function plotArduinoSig_journal()
	clc;
	close all;
	dir = "//DAXZ-SERVER/ResearchServer/ContactFreePulseMeasurement/Data/LAPSUN/S20/Video_Frame";
	plotInd = 1;
	timeRecord = 10;
	for file = {"50-normal", "50-press", "180-press"} 
		normal = load([dir "/" file{1} "_XoverY_ori/_rawPulseArduino_Sync.txt"]);
		time = normal(:, 1)*3600 + normal(:, 2)*60 + normal(:, 3);
		t0 = time(1) + 10
		t1 = t0+timeRecord
		ind0 = find(time==t0)(1);
		ind1 = find(time==t1)(1);
		normal = normal(ind0:ind1, 4);
		normal = [linspace(1, timeRecord, length(normal))', normal];
		subplot(3,1,plotInd); 
		hold on;
		plot(normal(:, 1), normal(:, 2), 'b')

		press = load([dir "/" file{1} "_XoverY_ori/_rawPulseArduino_Sync_press.txt"]);
		% time = press(:, 1)*3600 + press(:, 2)*60 + press(:, 3);
		% t0 = time(1000);
		% t1 = t0+timeRecord;
		ind0 = find(time==t0)(1);
		ind1 = find(time==t1)(1);
		press = press(ind0:ind1, 4);
		press = [linspace(1, timeRecord, length(press))', press];
		plot(press(:, 1), press(:, 2), 'r')
		xlim([1 timeRecord])
		hold off;
		if (plotInd == 1)
			h = legend('normal side', 'pressed side');
			rect = [0.75, 0.9, .15, .08];
			set(h, 'Position', rect)
		end
		if (plotInd == 2)
			ylabel('Amplitude')
		end
		plotInd = plotInd+1;
	end
	xlabel('Time(sec)')
	saveas(1, ["C:/Users/DAXZ-GAMING/Desktop/1/_arduino.png"]);
end