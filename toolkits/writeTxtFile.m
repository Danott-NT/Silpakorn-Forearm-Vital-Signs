function writeTxtFile(data_path, data)
	% data = rand(20000, 6);
	% data_path = "//DAXZ-SERVER/ResearchServer/ContactFreePulseMeasurement/Data/LAPSUN/S1/Feature/50-normal_EVM_Skin/data.txt"
	% tic
	fid = fopen(data_path, 'w+');
	for i=1:size(data, 1)
	    fprintf(fid, '%f\t', data(i,:));
	    fprintf(fid, '\n');
	end
	fclose(fid);
	% toc

end