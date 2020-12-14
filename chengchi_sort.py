f = h5py.File(file_to_load, 'r')
    num_channels = f['mat'].shape[0]
    ts2 = np.array(f['mat'])
    ts = np.transpose(ts2)
    f.close()
    sampling_frequency = 40000  # in Hz
    #timeseries = np.transpose(ts)
    geom = np.zeros((num_channels, 2))
    geom[:, 0] = range(num_channels)
    recording = se.NumpyRecordingExtractor(timeseries=ts, geom=geom, sampling_frequency=sampling_frequency)
10:15
# bandpass
    import spikeinterface.toolkit as st
    recording_f = st.preprocessing.bandpass_filter(recording, freq_min=300, freq_max=6000)
    import spikeinterface.sorters as ss
    default_ms4_params = ss.Mountainsort4Sorter.default_params()
    default_ms4_params['detect_threshold'] = 4
    default_ms4_params['curation'] = False
    sorting_MS4 = ss.run_mountainsort4(recording=recording_f, **default_ms4_params, output_folder='tmp_MS4/'+fileone)
    wf = st.postprocessing.get_unit_waveforms(recording, sorting_MS4, ms_before=1, ms_after=2,
                                              save_as_features=True, verbose=True)
    max_chan = st.postprocessing.get_unit_max_channels(recording, sorting_MS4, save_as_property=True, verbose=True)
    savename = '/Users/chengchichu/Desktop/ml4test/'+'bla_wvf_'+fileone.replace('.pl2','.mat')
    #scipy.io.savemat(savename,{'wf':wf,'maxchn':max_chan}, do_compression=True)
    #savename = 'test.mat'
    #a = {}
    #a['wf'] = wf
    hdf5storage.savemat(savename, {'wf':wf})