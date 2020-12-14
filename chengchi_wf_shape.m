list = dir('/Users/chengchichu/Desktop/ml4test/bla_wvf_*');
list2 = {list.name};
% mda path
pl2names = load('/Users/chengchichu/Desktop/dictator_analysis/pl2s_used.mat');
% load('toyset.mat')
% for iS = 4
cnt2 = 0;
pk_dist = [];
fr = [];
imaxs_ss = [];
tmplates = [];
for iS = 1:numel(list2)
  if iS == 45; continue, end
% for iS = 36
  list2{iS}
  mdapath = fullfile('/Users/chengchichu/tmp_MS4/', pl2names.session_used{iS},'firings.mda')
  spktime = readmda(mdapath);
  wvout = load(fullfile('/Users/chengchichu/Desktop/ml4test/',list2{iS}));
  try 
  assert(size(spktime,2) == sum(cellfun(@(x) size(x,1), wvout.wf)))
  catch err 
  assert(size(spktime,2) == size(wvout.wf,1))
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
  nUnits_recoveredid = unique(spktime(3,:));
  cnt2 = cnt2+nUnits_recovered;
  figure(1),clf,set(1,'Position',[37 156 998 42*(nCh)])
  cnt = 0;
  rmap = reshape(1:nCh*nUnits_recovered,[nUnits_recovered nCh])';
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
%        h = subplot(nCh,nUnits_recovered,rmap(cnt));
%        hs = [hs h];
       unit_wv = squeeze(wf{j}(:,i,:));
       template = mean(zscore(unit_wv,0,2));  
%       plot(template,'k','LineWidth',1)
%       axis off
       tmp_tmp{i} = template;
     end
     [~, i_varmax] = max(cellfun(@var,tmp_tmp));  
%     plot(tmp_tmp{i_varmax})
     [peak, i2] = min(tmp_tmp{i_varmax});  
     [trough, i1] = max(tmp_tmp{i_varmax}(i2:end));
     pkt = abs(i1+i2-1 - i2);
%     pkt = [i1 i2];
     mfr = length(find(spktime(3,:) == nUnits_recoveredid(j)))/ (spktime(2,end)/40000);
     pkts{j} = pkt;
     mfrs{j} = mfr;
     imaxs{j} = i_varmax;
     tmpss{j} = tmp_tmp;
%      linkaxes(hs,'xy');
  end
   pk_dist{iS} = pkts;
   fr{iS} = mfrs;
   imaxs_ss{iS} = imaxs;
   tmplates{iS} = tmpss;
  save_name = strrep(list2{iS},'.mat','');
%    printfigs('/eps',1,save_name,'png')
end