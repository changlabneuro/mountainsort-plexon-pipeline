function spike_compare = compare_spiketimes(old_spiketimes, ...
  new_spiketimes, sampling_rate)

% Function to compare information between the newly sorted units through
% mountainsort4 compared to the previously sorted units in mountainsort3
%
%
% Edit history:
% -------------
% % Feb 04, 2021 - Script made | Prabaha
% %

spike_compare = struct();

extracted_spike_info = extract_spike_info( old_spiketimes, new_spiketimes, sampling_rate );
unit_diff_per_file = get_unit_difference( extracted_spike_info );
spike_ct_diff_per_file = get_spike_ct_difference(extracted_spike_info);
isi_dist_t_test = compare_isi_dist_t_test( extracted_spike_info );
diff_between_spikes = find_diff_between_spikes(extracted_spike_info);

spike_compare.extracted_spike_info         = extracted_spike_info;
spike_compare.unit_diff_per_file           = unit_diff_per_file;
spike_compare.spike_ct_diff_per_file       = spike_ct_diff_per_file;
spike_compare.isi_dist_t_test              = isi_dist_t_test;
spike_compare.diff_between_spikes              = diff_between_spikes;

end

function extracted_spike_info = extract_spike_info( old_spiketimes, new_spiketimes, sampling_rate)

old_filename_list = {};
extracted_spike_info = {};



for old_file_ind = 1:numel( old_spiketimes.all_spike_time )
  old_filename_list{old_file_ind} = old_spiketimes.all_spike_time{old_file_ind}.filename{1}(1:end-4);
end

for new_file_ind = 1:numel( new_spiketimes.all_spiketimes )
  new_file_name = new_spiketimes.all_spiketimes{new_file_ind}{2};
  new_spikedata = new_spiketimes.all_spiketimes{new_file_ind}{1};
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
    new_units_data = extracted_spike_info{file_ind}.new_units_data.data;

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
    new_units_data = extracted_spike_info{file_ind}.new_units_data.data;
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

function new_data = get_adjusted_new_units_data(new_spikedata, sampling_rate)

if nargin < 2
  sampling_rate = 40000;
end
new_data = struct();
new_file_data = {};
spike_ct = [];
for unit = 1:numel(new_spikedata)
    spike_ct(unit) = numel(new_spikedata{unit});
end
sorted_spike_ct = sort(spike_ct, 'descend');
sorted_unit_id = {};
for unit = 1:numel(new_spikedata)
    
    idx = find(spike_ct == sorted_spike_ct(unit));
    sorted_unit_id{unit} = idx(1);
    new_file_data{unit} = new_spikedata{idx(1)} / sampling_rate;
end
new_data.unit_list = sorted_unit_id;
new_data.data = new_file_data;

end


function old_units_data = get_old_units_data(old_spikedata, old_file_ind)


old_units_data = {};

old_units_cellarray = old_spikedata{old_file_ind}.data;

spike_ct = [];
for unit = 1:numel(old_units_cellarray)
    spike_ct(unit) = numel(old_units_cellarray{unit}.data);
end
sorted_spike_ct = sort(spike_ct, 'descend');

for unit = 1:numel(old_units_cellarray)
    idx = find(spike_ct == sorted_spike_ct(unit));
    old_units_data{unit} = old_units_cellarray{idx}.data;
end

end

function diff_between_spikes = find_diff_between_spikes(extracted_spike_info)

diff_between_spikes = {};

for file_ind = 1:numel(extracted_spike_info)

    old_units_data = extracted_spike_info{file_ind}.old_units_data;
    new_units_data = extracted_spike_info{file_ind}.new_units_data.data;
    unit_comp = {};
    for new_unit_ind = 1:numel(new_units_data)
      for old_unit_ind = 1:numel(old_units_data)
          old_unit = old_units_data{old_unit_ind};
          new_unit = new_units_data{new_unit_ind};
          difference_closest = find_diff_with_closest_unit(old_unit, new_unit);
          unit_comp{old_unit_ind, new_unit_ind} = difference_closest;
      end
    end
    diff_between_spikes{file_ind} = unit_comp;
end

end



function difference_closest = find_diff_with_closest_unit(old_unit, new_unit)
  last_spike_ind = 1;
  difference_closest = [];
  start_new = new_unit(1);
  end_new = new_unit(end);
  start_old = old_unit(1);
  end_old = old_unit(end);
  if (start_new < start_old && end_new > start_old) ||...
          (start_new > start_old && end_old > start_new)
      for new_spike_ind = 1:numel(new_unit)
        min_diff = intmax;
        for old_spike_ind = last_spike_ind : numel(old_unit)
            diff_spikes = new_unit(new_spike_ind) - old_unit(old_spike_ind);
            if diff_spikes < min_diff
                min_diff = diff_spikes;
                last_spike_ind = old_spike_ind;
            else
                break   
            end 
        end
        difference_closest(new_spike_ind) = min_diff;
      end
  end
    

end