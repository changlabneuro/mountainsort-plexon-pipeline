#!/usr/bin/env python

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np 
import spikeinterface.comparison as sc
import spikeinterface.extractors as se
import spikeinterface.sorters as ss
import spikeinterface.toolkit as st
import spikeinterface.widgets as sw
import scipy.io
import h5py
import hdf5storage
import os
from sklearn import preprocessing 
from numpy.random import randint
from matplotlib.pyplot import cm


class Metrics():
  def __init__(self, file_name, metrics):
    self.file_name = file_name
    self.metrics = metrics

def extract_recording(timeseries):
  sampling_frequency = 40000
  geom = np.zeros((num_channels, 2))
  geom[:, 0] = range(num_channels)
  recording = se.NumpyRecordingExtractor(timeseries=timeseries, geom=geom, sampling_frequency=sampling_frequency)
  return recording

def preprocess_recording(recording, freq_min, freq_max):
  recording_f = st.preprocessing.bandpass_filter(recording, freq_min=freq_min, freq_max=freq_max)
  return recording_f


def sort_recording(recording, file_name):
  output_dir = '../tmp_MS4/'
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  default_ms4_params = ss.Mountainsort4Sorter.default_params()
  default_ms4_params['detect_threshold'] = 4
  default_ms4_params['curation'] = False
  default_ms4_params['filter'] = False
  sorting = ss.run_mountainsort4(recording=recording, **default_ms4_params, output_folder=output_dir+file_name)
  return sorting

def postprocess_recording(recording, sorting, file_name):
  output_dir = '../average_waveform/'
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)

  wf = st.postprocessing.get_unit_waveforms(recording, sorting, ms_before=1, ms_after=2,
                                            save_as_features=True, verbose=True)
  print("Shared unit spike feature names: ", sorting.get_shared_unit_spike_feature_names())
  print("Waveforms (n_spikes, n_channels, n_points)", wf[0].shape)

  max_chan = st.postprocessing.get_unit_max_channels(recording, sorting, save_as_property=True, verbose=False)
  savename = output_dir+file_name+'_avg_waveforms.mat'

  templates = st.postprocessing.get_unit_templates(recording, sorting, max_spikes_per_unit=200,
                                                  save_as_property=True, verbose=False)
    
  metrics= st.validation.compute_quality_metrics(sorting=sorting, recording=recording,
                                                metric_names=['firing_rate', 'isi_violation', 'snr', 'nn_hit_rate', 'nn_miss_rate'],
                                                as_dataframe=True)

  scipy.io.savemat(savename,{'wf':wf,'maxchn':max_chan, 'templates':templates, 'metrics':metrics}, do_compression=True)

  return wf, max_chan, templates, metrics

def visualize_units(recording, recording_f, sorting, wf, templates, num_channels, file_name):
    wf2 = st.postprocessing.get_unit_waveforms(recording, sorting, ms_before=4, ms_after=0,
                                          verbose=True)
    color= cm.rainbow(np.linspace(0,1,len(wf)))
    output_dir = "waveform_visualization/" + file_name +"/"
    if not os.path.exists(output_dir):
      os.mkdir(output_dir)
    fig, axs = plt.subplots(3, 2, constrained_layout = True)
    axe = axs.ravel()
    sw.plot_timeseries(recording_f, ax = axe[0])
    for i in range(len(wf)):
        normalized = preprocessing.normalize(np.array(templates)[i][:,:])
        axe[1].plot(normalized.T, label=i+1, color=color[i])
    sw.plot_rasters(sorting, ax = axe[2])
    sw.plot_pca_features(recording, sorting, colormap='rainbow', nproj=3, max_spikes_per_unit=100, axes = axe[3:6])
    legend = axe[1].legend(loc='lower right')
    fig.set_size_inches((8.5, 11), forward=False)
    # name = file_name.split("_")
    # # acc_1_04052016_kurocoppola_pre
    fig.suptitle(file_name)
    fig.savefig(output_dir +'_full_info.png', dpi=500)
    for idx in range(len(wf)):
      for j in range(num_channels):
        fig, axs = plt.subplots(4, 3, constrained_layout = True)
        axe = axs.ravel()
        normalized = preprocessing.normalize(wf[idx][:,j,:])
        axe[2].plot(normalized.T, lw=0.3, alpha=0.3)
        normalized = preprocessing.normalize(np.array(templates)[idx][:,:])
        axe[2].plot(normalized.T, lw=0.8, color='black')
        axe[0].plot(templates[idx].T, lw=0.8)
        for i in range(len(wf)):
            normalized = preprocessing.normalize(np.array(templates)[i][:,:])
            if idx == i:
                axe[1].plot(normalized.T, label=i+1)
            else:
                axe[1].plot(normalized.T, label=i+1, alpha=0.3)
        values = randint(0, len(wf), 3)
        for i, val in enumerate(values):    
            axe[3 + i].plot(wf2[idx][val, j, :].T)
        sw.plot_pca_features(recording, sorting, colormap='binary', figure = fig, nproj=3, axes = axe[6:9], max_spikes_per_unit=100)
        sw.plot_pca_features(recording, sorting, unit_ids = [idx], figure = fig, nproj=3, axes = axe[6:9] ,max_spikes_per_unit=100)
        sw.plot_isi_distribution(sorting, unit_ids = [idx + 1], bins=100, window=1, axes = axe[9: 11])
        sw.plot_amplitudes_distribution(recording, sorting, max_spikes_per_unit=300, unit_ids = [idx + 1], axes = axe[10: 12])
        sp = sorting.get_unit_spike_train(unit_id = idx + 1)
        axe[11].plot(sp)
        legend = axe[1].legend(loc='lower right', shadow=True)
        axe[0].set_title('Template')
        axe[2].set_title('Waveforms, template')
        axe[3].set_title('Random waveforms')
        axe[6].set_title('PCA component analysis')
        axe[9].set_title('ISI distribution')
        axe[10].set_title('Amplitude distribution')
        axe[11].set_title('Spike frame')
        fig.set_size_inches((8.5, 11), forward=False)
        fig.suptitle("File name: " + file_name + '\nChannel: ' + str(j+1).zfill(2) + '\nUnit: ' + str(idx + 1).zfill(2))
        fig.savefig(output_dir + file_name + '_ch' + str(j+1).zfill(2) + '_u'+ str(idx+1).zfill(2) + '.png', dpi=500)
    plt.close('all')

if __name__ == "__main__":
  # directory_in_str = input("Input directory location here")
  directory_in_str = "../eisi_raw_mat_files"
  directory = os.fsencode(directory_in_str)
  sessions_used = []
  metrics_sessions = []
  for file in os.listdir(directory):
    file_name = os.fsdecode(file)
    # file_name = "acc_1_04052016_kurocoppola_pre.mat"
    file_path = os.path.join(directory_in_str, file_name)
    file_name = file_name.replace('.mat', '')
    print(file_name)
    f = h5py.File(file_path, 'r')
    num_channels = f['mat'].shape[1]
    ts = np.transpose(np.array(f['mat']))
    f.close()
    recording = extract_recording(timeseries=ts)
    recording_f = preprocess_recording(recording=recording, freq_min=300, freq_max=6000)
    sorting = sort_recording(recording=recording_f, file_name=file_name)
    wf, max_chan, templates, metrics = postprocess_recording(recording = recording, sorting = sorting, file_name = file_name)
    output_dir = "waveform_visualization/"
    if not os.path.exists(output_dir):
      os.mkdir(output_dir)
    visualize_units(recording, recording_f, sorting, wf, templates, num_channels, file_name)
    if metrics.shape[1] == 5:
      metrics_sessions.append(Metrics(file_name, metrics))
    sessions_used.append(file_name)
  scipy.io.savemat('metrics.mat',{'sessions_used':sessions_used, 'metrics_sessions': metrics_sessions}, do_compression=True)
    

