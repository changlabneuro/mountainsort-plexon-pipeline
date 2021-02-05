function extracted_spike_info = compare_spiketimes(old_spiketimes, ...
  new_spiketimes)

% Function to compare information between the newly sorted units through
% mountainsort4 compared to the previously sorted units in mountainsort3
%
%
% Edit history:
% -------------
% % Feb 04, 2021 - Script made | Prabaha
% %

extracted_spike_info = extract_spike_info( old_spiketimes, new_spiketimes );
unit_diff_pre_file = get_unit_difference( extracted_spike_info );

end

function extracted_spike_info = extract_spike_info( old_spiketimes, new_spiketimes )

extracted_spike_info = {};

for new_file_ind = 1:numel( new_spiketimes )
  new_file_name = new_spiketimes.all_spiketimes{new_file_ind}{2};
  new_file_data = new_spiketimes.all_spiketimes{new_file_ind}{1};
  old_spikedata = old_spiketimes.all_spike_time;
  old_file_ind = find( strcmp( new_file_name, old_spikedata.filename(1:end-4) ) );
  if isempty(old_file_ind)
    fprintf('Error! No old file by the name : %s\n', new_file_name);
  elseif numel(old_file_ind > 1)
    fprintf('Error! More than one old file named : %s\n', new_file_name);
  else
    extracted_spike_info{end+1} = make_spike_info();
    old_file_data = get_old_units_data( old_spikedata, old_file_ind );
    extracted_spike_info.old_units_data = old_file_data;
    extracted_spike_info.num_old_units = numel( old_file_data );
    extracted_spike_info.new_units_data = new_file_data;
    extracted_spike_info.num_new_units = numel( new_file_data );
  end
end

end

function unit_diff_pre_file =  get_unit_difference(extracted_spike_info)

num_units_per_file_new = [ extracted_spike_info{1:end}.num_old_units ];
num_units_per_file_new = [ extracted_spike_info{1:end}.num_new_units ];

unit_diff_pre_file = num_units_per_file_new - num_units_per_file_new;

end

function extracted_spike_info = make_spike_info()

extracted_spike_info = struct();
extracted_spike_info.old_units_data = nan;
extracted_spike_info.num_old_units = nan;
extracted_spike_info.new_units_data = nan;
extracted_spike_info.num_new_units = nan;

end


function old_units_data = get_old_units_data(old_spikedata, old_file_ind)

old_units_data = old_spikedata{old_file_ind}.data;

end