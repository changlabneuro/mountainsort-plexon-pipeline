function spike_compare = compare_spiketimes(old_spiketimes, ...
  new_spiketimes)

% Function to compare information between the newly sorted units through
% mountainsort4 compared to the previously sorted units in mountainsort3
%
%
% Edit history:
% -------------
% % Feb 04, 2021 - Script made | Prabaha
% %

spike_compare = struct();

extracted_spike_info = extract_spike_info( old_spiketimes, new_spiketimes );
unit_diff_per_file = get_unit_difference( extracted_spike_info );
spike_ct_diff_per_file = get_spike_ct_difference(extracted_spike_info);
isi_dist_t_test = compare_isi_dist_t_test( extracted_spike_info );

spike_compare.extracted_spike_info         = extracted_spike_info;
spike_compare.unit_diff_per_file           = unit_diff_per_file;
spike_compare.spike_ct_diff_per_file       = spike_ct_diff_per_file;
spike_compare.isi_dist_t_test              = isi_dist_t_test;

end

function extracted_spike_info = extract_spike_info( old_spiketimes, new_spiketimes )

old_filename_list = {};
extracted_spike_info = {};

for old_file_ind = 1:numel( old_spiketimes.all_spike_time )
  old_filename_list{old_file_ind} = old_spiketimes.all_spike_time{old_file_ind}.filename{1}(1:end-4);
end

for new_file_ind = 1:numel( new_spiketimes.all_spiketimes )
  new_file_name = new_spiketimes.all_spiketimes{new_file_ind}{2};
  new_spikedata = new_spiketimes.all_spiketimes{new_file_ind}{1};
  sampling_rate = 40000;
  %   new_file_data = new_spikedata if using original sample number
  new_file_data = get_adjusted_new_units_data(new_spikedata, sampling_rate);
  
  old_spikedata = old_spiketimes.all_spike_time;
  old_file_ind = find( strcmp( new_file_name, old_filename_list ) );
  if isempty(old_file_ind)
    fprintf('Error! No old file by the name : %s\n', new_file_name);
  elseif numel(old_file_ind) > 1
    fprintf('Error! More than one old file named : %s\n', new_file_name);
  else
    extracted_spike_info{end+1} = make_spike_info();
    old_file_data = get_old_units_data( old_spikedata, old_file_ind );
    extracted_spike_info{end}.filename = new_file_name;
    extracted_spike_info{end}.old_units_data = old_file_data;
    extracted_spike_info{end}.num_old_units = numel( old_file_data );
    extracted_spike_info{end}.new_units_data = new_file_data;
    extracted_spike_info{end}.num_new_units = numel( new_file_data );
    
  end
end

end

function unit_diff_per_file =  get_unit_difference(extracted_spike_info)

spike_info_array = [  extracted_spike_info{1:end} ];


num_units_per_file_old = [ spike_info_array.num_old_units ];
num_units_per_file_new = [ spike_info_array.num_new_units ];

unit_diff_per_file = num_units_per_file_new - num_units_per_file_old;

end

function spike_ct_diff_per_file = get_spike_ct_difference(extracted_spike_info)

spike_ct_diff_per_file = {};

for file_ind = 1:numel(extracted_spike_info)
    
    file_ct_diff = {};

    old_units_data = extracted_spike_info{file_ind}.old_units_data;
    new_units_data = extracted_spike_info{file_ind}.new_units_data;

    for new_unit_ind = 1:numel(new_units_data)
        
      new_unit_ct_diff = {};
      for old_unit_ind = 1:numel(old_units_data)
          old_unit_spike_ct = numel(old_units_data{old_unit_ind});
          new_unit_spike_ct = numel(new_units_data{new_unit_ind});
          new_unit_ct_diff{old_unit_ind} = new_unit_spike_ct - old_unit_spike_ct;
      end
      file_ct_diff{new_unit_ind} = new_unit_ct_diff;
    end
    spike_ct_diff_per_file{file_ind} = file_ct_diff;
end


end

function isi_dist_t_test = compare_isi_dist_t_test(extracted_spike_info)

isi_dist_t_test = {};

for file_ind = 1:numel(extracted_spike_info)
    
    file_z_test = {};

    old_units_data = extracted_spike_info{file_ind}.old_units_data;
    new_units_data = extracted_spike_info{file_ind}.new_units_data;
    unit_comp_p = [];
    unit_comp_h = [];
    for new_unit_ind = 1:numel(new_units_data)
      for old_unit_ind = 1:numel(old_units_data)
          old_unit_isi_dist = diff(old_units_data{old_unit_ind});
          new_unit_isi_dist = diff(new_units_data{new_unit_ind});
          [h, p] = ttest2(old_unit_isi_dist, new_unit_isi_dist);
          unit_comp_p(old_unit_ind, new_unit_ind) = p;
          unit_comp_h(old_unit_ind, new_unit_ind) = h;
      end
    end
    isi_dist_t_test{file_ind} = {unit_comp_h, unit_comp_p};
end

end

function extracted_spike_info = make_spike_info()

extracted_spike_info = struct();
extracted_spike_info.filename = [];
extracted_spike_info.old_units_data = nan;
extracted_spike_info.num_old_units = nan;
extracted_spike_info.new_units_data = nan;
extracted_spike_info.num_new_units = nan;

end

function new_file_data = get_adjusted_new_units_data(new_spikedata, sampling_rate)

if nargin < 2
  sampling_rate = 40000;
end
new_file_data = {};
for unit = 1:numel(new_spikedata)
  new_file_data{unit} = new_spikedata{unit}/sampling_rate;
end

end



function old_units_data = get_old_units_data(old_spikedata, old_file_ind)
old_units_data = {};

old_units_cellarray = old_spikedata{old_file_ind}.data;
for unit = 1:numel(old_units_cellarray)
  old_units_data{unit} = old_units_cellarray{unit}.data;
end

end