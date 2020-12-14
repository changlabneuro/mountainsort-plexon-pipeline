# Mountain Sort Documentation

## Extraction

- For MDA format, MountiainSort requires raw.mda file with M x N timeseries in mda format where M is number of channels and N is number of timepoints. In our case, we are using .mat files with same dimensions.

- It would also require geom.csv that represents the probe geometry and is a comma-separated text file containing 2-d or 3-d coordinates of the electrodes. The number of lines (or rows) in geom.csv should equal M, the number of channels. We are using geom[:, 0] = range(num_channels) for specifying the geometry.

- params.json is a JSON format file containing dataset-specific parameters, including samplerate (in Hz) and spike_sign (-1, 0, or 1).
sampling_frequency 40000 Hz

## Sorting
- default parameters: {'detect_sign': -1, 'adjacency_radius': -1, 'freq_min': 300, 'freq_max': 6000, 'filter': True, 'whiten': True, 'curation': False, 'num_workers': None, 'clip_size': 50, 'detect_threshold': 3, 'detect_interval': 10, 'noise_overlap_threshold': 0.15} 

- curation = false as small number of channels
- filter = false as using preprocessing filter
- detect_threshold = 4
- detect_sign = -1 as looking for negative spike peaks with most electrophysiology datasets are made up almost entirely of negative spikes
- adjacency_radius = -1 as using only one electrode neighborhood containing all the channels

## Post-processing

https://spikeinterface.readthedocs.io/en/latest/modules/toolkit/plot_2_postprocessing.html#sphx-glr-modules-toolkit-plot-2-postprocessing-py

- Assuming the sorting is the output of a spike sorter, the postprocessing module allows to extract all relevant information from the paired recording-sorting.

- Waveforms are extracted with the get_unit_waveforms function by extracting snippets of the recordings when spikes are detected. When waveforms are extracted, the can be loaded in the SortingExtractor object as features. The ms before and after the spike event can be chosen. Waveforms are returned as a list of np.arrays (n_spikes, n_channels, n_points)

- Similarly to waveforms, templates - average waveforms - can be easily  extracted using the get_unit_templates. When spike trains have  numerous spikes, you can set the max_spikes_per_unit to be extracted.  If waveforms have already been computed and stored as features, those  will be used. Templates can be saved as unit properties.

- In the same way, one can get the recording channel with the maximum  amplitude and save it as a property.

- For some applications, for example validating the spike sorting output,  PCA scores can be computed. PCA scores can be also computed electrode-wise. In the previous example,  PCA was applied to the concatenation of the waveforms over channels.
