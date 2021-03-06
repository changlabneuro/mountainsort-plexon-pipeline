function generate_spike_comparison_plots(old_spiketimes, new_spiketimes)

% Function to generate plots comparing information between the newly sorted 
% units through mountainsort4 compared to the previously sorted units in 
% mountainsort3
%
%
% Edit history:
% -------------
% % Feb 17, 2021 - Script made (just t-test) | Prabaha
% %

% if compare spiketime output is not saved

if nargin < 2
  old_spiketimes = load('dictator_game_SUAdata_pre.mat');
  new_spiketimes = load('ml4alg_brains_spiketimes_2021.mat');
end

spike_compare = util.compare_spiketimes(old_spiketimes, new_spiketimes, 40000);

% if compare spiketime output is saved as spiketimes_compare.mat file
% spike_compare = load('spiketimes_compare.mat');
% spike_compare = spike_compare.temp;

% plot_isi_dist_comp(spike_compare);

plot_spike_time_difference_histogram(spike_compare);

% plot_raster_plot(spike_compare);

end

function plot_raster_plot(spike_compare)

% if ~exist('/spike_comp_plots/raster_plots', 'dir')
%    mkdir('/spike_comp_plots/raster_plots/');
% end

save_path = './spike_comp_plots/raster_plots/';
spike_info = spike_compare.extracted_spike_info;

for file_ind = 1:numel(spike_info)
  old_units_data = spike_compare.extracted_spike_info{file_ind}.old_units_data;
  new_units_data = spike_compare.extracted_spike_info{file_ind}.new_units_data.data;
  
  figure();
  
  min_len = 500;
  for i = 1:numel(old_units_data)
      if size(old_units_data{i}, 2) < min_len; min_len = size(old_units_data{i}, 2); end
  end
  for i = 1:numel(new_units_data)
      if size(new_units_data{i}, 2) < min_len; min_len = size(new_units_data{i}, 2); end
  end
  old_new_spiketimes = {};
  trials = {};
  for i = 1:numel(old_units_data)
      curr_spike_time = old_units_data{i}(1:min_len);
      old_new_spiketimes{i} = curr_spike_time(1:min_len);
      trials{i}=i*ones(size(old_new_spiketimes{i}));
  end
  
  for i = numel(old_units_data) + 1 :numel(old_units_data) + numel(new_units_data)
    ind = i - numel(old_units_data);
    curr_spike_time = new_units_data{ind}(1:min_len);
    old_new_spiketimes{i} = curr_spike_time(1:min_len);
    trials{i}=i*ones(size(old_new_spiketimes{i}));
  end
  util.rasterplot(cell2mat(old_new_spiketimes(:)), cell2mat(trials(:)), numel(old_new_spiketimes));
  xlabel('Time (s)');
  ylabel(['Units, Old: 1 - ' num2str(i - numel(new_units_data)) ', New: ' num2str(i - numel(new_units_data)) ' - ' num2str(i)]);
  
  plot_filename = [save_path 'rasterplot' spike_info{file_ind}.filename '.pdf'];
  saveas(gcf, plot_filename);
  close all;
end


end

function plot_spike_time_difference_histogram(spike_compare)

% if ~exist('/spike_comp_plots/spike_time_diff_comp_plots', 'dir')
%    mkdir('/spike_comp_plots/spike_time_diff_comp_plots/');
% end

save_path = './spike_comp_plots/spike_time_diff_comp_plots/';
diff_between_spikes = spike_compare.diff_between_spikes;
spike_info = spike_compare.extracted_spike_info;

for file_ind = 1:numel(spike_info)
  spike_ct_diff = diff_between_spikes{file_ind};
  
  figure();
  
  num_old_units = size(spike_ct_diff, 1);
  num_new_units = size(spike_ct_diff, 2);
  
  for old_ind = 1:num_old_units
      for new_ind = 1:num_new_units
          subplot(num_old_units, num_new_units, old_ind + (new_ind-1)*num_old_units);
          histogram(spike_ct_diff{old_ind,new_ind});
          % new_unit_ind - unit index from sorted file, where new_ind is just index when sorted in ascending order      
          new_unit_ind = spike_info{file_ind}.new_units_data.unit_list{new_ind};
          title(['O: ' num2str(old_ind) ' N: ' num2str(new_ind) '-' num2str(new_unit_ind)]);
          
      end
  end
  sgtitle(['Spike Difference Between New and Old Units | ' spike_info{file_ind}.filename]);
  
  plot_filename = [save_path 'spike_diff' spike_info{file_ind}.filename '.pdf'];
  saveas(gcf, plot_filename);
  close;
  
end

end

function plot_isi_dist_comp(spike_compare)

if ~exist('spike_comp_plots', 'dir')
   mkdir('spike_comp_plots');
end
save_path = './spike_comp_plots/isi_dist_t_test/';
isi_dist_t_test = spike_compare.isi_dist_t_test;
spike_info = spike_compare.extracted_spike_info;


for file_ind = 1:numel(isi_dist_t_test)
  h_mat = isi_dist_t_test{file_ind}{1};
  p_mat = isi_dist_t_test{file_ind}{2};
  
  figure();
  subplot(1,2,1);
  imagesc(h_mat);
  colormap gray;
  colorbar;
  caxis([0 1]);
  title('t-test result');
  xlabel('New Unit Number');
  ylabel('Old Unit Number');
  
  subplot(1,2,2);
  imagesc(p_mat);
  colormap autumn;
  colorbar;
  for old_unit_ind = 1:size( p_mat, 1 )
    for new_ind = 1:size( p_mat, 2 )
      new_unit_ind = spike_info{file_ind}.new_units_data.unit_list{new_ind};
      text(new_unit_ind, old_unit_ind, num2str(p_mat(old_unit_ind,new_unit_ind),2), ...
        'FontSize', 7);
    end
  end
  title('t-test p-value');
  xlabel('New Unit Number');
  ylabel('Old Unit Number');
  
  suptitle(['ISI dist comp | t-test | ' spike_info{file_ind}.filename]);
  
  plot_filename = [save_path 'isi_dist_t_test-' spike_info{file_ind}.filename '.pdf'];
  saveas(gcf, plot_filename);
  close;
  
end

end