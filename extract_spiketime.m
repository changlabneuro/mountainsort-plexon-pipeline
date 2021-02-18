fid = fopen('pl2_used.txt','rt');
all_spiketimes = {};
idx = 1;
while true
  thisline = fgetl(fid);
  if ~ischar(thisline); break; end  %end of file
  mdapath = fullfile('/Users/beam/Documents/chang_lab/tmp_MS4/', thisline,'firings.mda');
  spktime = readmda(mdapath);
  unit_spiketimes = {};
  spiketimes = {};
  for i = unique(spktime(3,:))
     spktime_filtered  = spktime(:, spktime(3,:) == i);
     unit_spktime = spktime_filtered(2, :);
     unit_spiketimes{i} = unit_spktime;
     spiketimes{1} = unit_spiketimes;
     spiketimes{2} = thisline;  
  end
  all_spiketimes{idx} = spiketimes;
  idx = idx + 1;
end
fclose('pl2_used.txt');