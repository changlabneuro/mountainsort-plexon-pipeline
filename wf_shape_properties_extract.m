list = {dir('/Users/beam/Documents/chang_lab/average_waveform/acc_*').name};
% mda path
pl2names = load('/Users/beam/Documents/chang_lab/main/pl2_used.mat');
pl2_ct = 0;
peak_distance = [];
frequency = [];
imaxs_ss = [];
templates = [];
for iS = 2:2
  mdapath = fullfile('/Users/beam/Documents/chang_lab/tmp_MS4/', strrep(list{iS},'_avg_waveforms.mat',''),'firings.mda');
  spktime = readmda(mdapath);
%   mdaname = strtrim(pl2names.sessions_used(iS, :));
%   matname = list{iS};
  list{iS}
  wvout = load(fullfile('/Users/beam/Documents/chang_lab/average_waveform/',list{iS}));
  try 
  assert(size(spktime,2) == sum(cellfun(@(x) size(x,1), wvout.wf)))
  catch err 
  end
  if isempty(wvout.wf) == 1
  continue
  end
  if iscell(wvout.wf)==0
    wf = {wvout.wf};
  else
    wf = wvout.wf;
  end
  nCh = size(wf{1},2);
  nUnits_recovered = numel(wf);
  nUnits_recovered_id = unique(spktime(3,:));
  pl2_ct = pl2_ct + nUnits_recovered;
  figure(1),clf,set(1,'Position',[37 156 998 42*(nCh)])
  cnt = 0;
  rmap = reshape(1:nCh * nUnits_recovered,[nUnits_recovered nCh])';
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
       h = subplot(nCh,nUnits_recovered,rmap(cnt));
       hs = [hs h];
       unit_wv = squeeze(wf{j}(:,i,:));
       template = mean(zscore(unit_wv,0,2));  
       plot(template,'LineWidth',1)
       axis off
       tmp_tmp{i} = template;
     end
     [~, i_varmax] = max(cellfun(@var,tmp_tmp));  
     plot(tmp_tmp{i_varmax})
     [peak, i2] = min(tmp_tmp{i_varmax});  
     [trough, i1] = max(tmp_tmp{i_varmax}(i2:end));
     pkt = abs(i1+i2-1 - i2);
    pkt = [i1 i2];
     mfr = length(find(spktime(3,:) == nUnits_recovered_id(j)))/ (spktime(2,end)/40000);
     pkts{j} = pkt;
     mfrs{j} = mfr;
     imaxs{j} = i_varmax;
     tmpss{j} = tmp_tmp;
     linkaxes(hs,'xy');
  end
   peak_distance{iS} = pkts;
   frequency{iS} = mfrs;
   imaxs_ss{iS} = imaxs;
   templates{iS} = tmpss;
   save_name = strrep(list{iS},'.mat','');
   print(save_name,'-dpng');
end