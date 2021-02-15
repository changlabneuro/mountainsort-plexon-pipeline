fid = fopen('pl2_used.txt','rt');
idx = 1;
peak_distance = [];
frequency = [];
imaxs_ss = [];
templates = [];
file_names = [];
while true
  thisline = fgetl(fid);
  if ~ischar(thisline); break; end  %end of file
  mdapath = fullfile('/Users/beam/Documents/chang_lab/tmp_MS4/', thisline,'firings.mda');
  wvout = load(strcat('/Users/beam/Documents/chang_lab/average_waveform/',thisline, '_avg_waveforms.mat'));
  spktime = readmda(mdapath);
  if iscell(wvout.wf)==0
    wf = {wvout.wf};
  else
    wf = wvout.wf;
  end
  nCh = size(wf{1},2);
  nUnits_recovered = numel(wvout.maxchn);
  nUnits_recovered_id = unique(spktime(3,:));
  cnt = 0;
  pkts = {};
  mfrs = {};
  imaxs = {};
  tmps = {};
  for j = 1:(nUnits_recovered)
     tmp_tmp = wvout.templates(j, wvout.maxchn(j) + 1, :);
     [peak, i2] = max(tmp_tmp);  
     [trough, i1] = min(tmp_tmp(i2:end));
     pkt = [peak trough peak-trough i2 i1 abs(i2 - i1)];
     mfr = length(find(spktime(3,:) == nUnits_recovered_id(j)))/ (spktime(2,end)/40000);
     pkts{j} = pkt;
     mfrs{j} = mfr;
     imaxs{j} = i_varmax;
     tmps{j} = tmp_tmp;
   end
   peak_distance{idx} = pkts;
   frequency{idx} = mfrs;
   imaxs_ss{idx} = imaxs;
   templates{idx} = tmps;
   file_names{idx} = thisline;
   idx = idx + 1;
end
fclose('all');