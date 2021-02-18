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

if nargin < 1
  old_spiketimes = load('dictator_game_SUAdata_pre.mat');
  new_spiketimes = load('ml4alg_brains_spiketimes_2021.mat');
end

spike_compare = util.compare_spiketimes(old_spiketimes, new_spiketimes);

plot_isi_dist_comp(spike_compare);

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
    for new_unit_ind = 1:size( p_mat, 2 )
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