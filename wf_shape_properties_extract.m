fid = fopen('pl2_used.txt','rt');
idx = 1;
peak_distance = [];
frequency = [];
imaxs_ss = [];
templates = [];
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
  nUnits_recovered = numel(wf);
  nUnits_recovered_id = unique(spktime(3,:));
  cnt = 0;
  pkts = {};
  mfrs = {};
  imaxs = {};
  tmpss = {};
  for j = 1:(nUnits_recovered)  
     hs = [];
     tmp_tmp = {};
     for i = 1:nCh
       idx_real = find(wvout.maxchn==i-1);  
       cnt = cnt+1;  
       unit_wv = squeeze(wf{j}(:,i,:));
       template = mean(zscore(unit_wv,0,2));  
       tmp_tmp{i} = template;
     end
     [~, i_varmax] = max(cellfun(@var,tmp_tmp));  
     [peak, i2] = min(tmp_tmp{i_varmax});  
     [trough, i1] = max(tmp_tmp{i_varmax}(i2:end));
     pkt = [abs(i2 - i1) i2 i1];
     mfr = length(find(spktime(3,:) == nUnits_recovered_id(j)))/ (spktime(2,end)/40000);
     pkts{j} = pkt;
     mfrs{j} = mfr;
     imaxs{j} = i_varmax;
     tmpss{j} = tmp_tmp;
   end
   peak_distance{idx} = pkts;
   frequency{idx} = mfrs;
   imaxs_ss{idx} = imaxs;
   templates{idx} = tmpss;
   idx = idx + 1;
end
fclose('pl2_used.txt');