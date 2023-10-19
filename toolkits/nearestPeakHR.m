function [hr, mag] = nearestPeakHR(ySpecCell, specIndCell, peakIndCell, pulseRefCell)
  ySpecCell = cellfun(@(ySpec, peakInd) arrayfun(@(channel) ySpec(:, channel)(peakInd{channel})', 1:size(ySpec, 2), 'UniformOutput',false), ySpecCell, peakIndCell, 'UniformOutput',false);
  ySpecCell = cellfun(@(ySpec) cell2mat(ySpec), ySpecCell, 'UniformOutput',false);
  specIndCell = cellfun(@(specInd, peakInd) specInd(cell2mat(peakInd))', specIndCell, peakIndCell, 'UniformOutput',false);

  [m, ind] = arrayfun(@(windowInd) min(abs(specIndCell{windowInd}-pulseRefCell(windowInd))), 1:size(pulseRefCell, 2));
  
  hr = arrayfun(@(windowInd) specIndCell{windowInd}(ind(windowInd)), 1:size(ind, 2));
  mag = arrayfun(@(windowInd) ySpecCell{windowInd}(ind(windowInd)), 1:size(ind, 2));
end