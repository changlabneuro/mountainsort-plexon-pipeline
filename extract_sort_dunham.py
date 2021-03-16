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
import logging
import pandas as pd
import math
from matplotlib.lines import Line2D

def extract_recording(timeseries, num_channels):
  sampling_frequency = 40000 # in Hz
  geom = np.zeros((num_channels, 2)) # zeros for channel relative location, input actual coordinates if known
  geom[:, 0] = range(num_channels)
  recording = se.NumpyRecordingExtractor(timeseries=timeseries, geom=geom, sampling_frequency=sampling_frequency)
  return recording

def preprocess_recording(recording, freq_min, freq_max):
  # apply bandpass filter with min and max frequency
  recording_f = st.preprocessing.bandpass_filter(recording, freq_min=freq_min, freq_max=freq_max)
  return recording_f


def sort_recording(recording, file_name):
  # modify default params, turn filter off if already filtered
  output_dir = '../tmp_MS4/'
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  default_ms4_params = ss.Mountainsort4Sorter.default_params()
  default_ms4_params['detect_threshold'] = 4
  default_ms4_params['curation'] = False
  default_ms4_params['filter'] = True
  sorting = ss.run_mountainsort4(recording=recording, **default_ms4_params, output_folder=output_dir+file_name)
  return sorting

def postprocess_recording(recording, sorting, file_name):
  # extract waveforms for units sorted with spike maximum channel, templates/ average waveforms and metrics values
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
  # uses firing rate, isi violation, snr for metrics values computed, can add more like nn_hit_rate and miss_rate
  metrics= st.validation.compute_quality_metrics(sorting=sorting, recording=recording,
                                                metric_names=['firing_rate', 'isi_violation', 'snr'],
                                                as_dataframe=True)

  scipy.io.savemat(savename,{'wf':wf,'maxchn':max_chan, 'templates':templates, 'metrics':metrics}, do_compression=True)

  return wf, max_chan, templates, metrics

def visualize_units(recording, recording_f, sorting, wf, templates, max_chan, num_channels, file_name):
  # Create color map for legend
  color= cm.rainbow(np.linspace(0,1,len(wf)))
  custom_lines = []
  for i in range(len(wf)):
    custom_lines.append(Line2D([0], [0], color=color[i], lw=1))
  # Create directory for file under waveform_visualization
  output_dir = "waveform_visualization/" + file_name +"/"
  if not os.path.exists(output_dir):
    os.mkdir(output_dir)
  # Plot general full information including timeseries data of bandpass filtered recording, templates of units with channel showing highest amplitude difference
  # raster information and 3 component PCA score analysis
  fig, axs = plt.subplots(3, 2, constrained_layout = True)
  axe = axs.ravel()
  sw.plot_timeseries(recording_f, ax = axe[0])
  # Normalize templates and plot for channel with max difference
  normalized_templates = []
  for i in range(len(wf)):     
      normalized_templates.append(preprocessing.normalize(np.array(templates)[i][:,:]))
      axe[1].plot(normalized_templates[i][max_chan[i]].T, color=color[i], label= i + 2)
  sw.plot_rasters(sorting, ax = axe[2])
  

  # Calculate pca scores with 3 component analysis and plot for multichannel recordings or use pre built one for single channel ones
  # Comment this section out and adjust the figure axes if taking to long to compute or not needed
  


  '''
  if num_channels > 1:
    pca_scores = st.postprocessing.compute_unit_pca_scores(recording, sorting, unit_ids = list(range(1, len(wf) + 1)), channel_ids = max_chan, n_comp=3, verbose=True)
    for i in range(len(wf)):
        axe[3].plot(pca_scores[i][:, max_chan[i], 0], pca_scores[i][:, max_chan[i], 1], '*')
        axe[4].plot(pca_scores[i][:, max_chan[i], 0], pca_scores[i][:, max_chan[i], 2], '*')
        axe[5].plot(pca_scores[i][:, max_chan[i], 1], pca_scores[i][:, max_chan[i], 2], '*')
  else:
    sw.plot_pca_features(recording, sorting, colormap='rainbow', nproj=3, max_spikes_per_unit=100, axes = axe[3:6])
  axe[1].legend(custom_lines, range(len(wf)), ncol=2)
  fig.set_size_inches((8.5, 11), forward=False)
  fig.suptitle(file_name + '\nNum of channels: ' + str(num_channels) + 'Num of units: ' + str(len(wf)) +'\nMax. channels: ' + str([x + 1 for x in max_chan]) + '\nMetrics valid units: ' + str(metrics.index[metrics['isi_violation'] < 1.5].tolist()))
  fig.savefig(output_dir +'_full_info.png', dpi=500)
  plt.close()

  # Plot pca score analysis for units with isi violations less than 1.5 with pca analysis using built in function and plotting pca scores after analysis for multichannel ones
  if num_channels > 1:
      fig, axs = plt.subplots(2, 3, constrained_layout = True)
  else:
      fig, axs = plt.subplots(1, 3, constrained_layout = True)
  sw.plot_pca_features(recording, sorting, unit_ids = metrics.index[metrics['isi_violation'] < 1.5].tolist(), figure = fig, nproj=3, axes = axs[0:3] ,max_spikes_per_unit=100)
  if num_channels > 1:
      for i in metrics.index[metrics['isi_violation'] < 1.5].tolist():
          axs[3].plot(pca_scores[i-1][:, max_chan[i-1], 0], pca_scores[i-1][:, max_chan[i-1], 1], '*')
          axs[4].plot(pca_scores[i-1][:, max_chan[i-1], 0], pca_scores[i-1][:, max_chan[i-1], 2], '*')
          axs[5].plot(pca_scores[i-1][:, max_chan[i-1], 1], pca_scores[i-1][:, max_chan[i-1], 2], '*', label=i+1)
      axe[5].legend(custom_lines, range(len(wf)), ncol=2)
  fig.set_size_inches((8.5, 11), forward=False)
  fig.suptitle(file_name + '\nNum of channels: ' + str(num_channels) + ' Num of units: ' + str(len(wf)) +'\nMax. channels: ' + str([max_chan[i-1] for i in metrics.index[metrics['isi_violation'] < 1.5].tolist()]) + '\nMetrics valid units: ' + str(metrics.index[metrics['isi_violation'] < 1.5].tolist()))
  fig.savefig(output_dir + file_name + '_pca_analysis.png', dpi=500)
  plt.close()

  # Plot template, template among other unit templates, template among waveform visualization for the unit, 3 randomly selected waveform visualizations, PCA score analysis with unit highlighted
  # isi and amplitude distribution and finally spike train per unit
  for idx in range(len(wf)):
    fig, axs = plt.subplots(4, 3, constrained_layout = True)
    axe = axs.ravel()
    # plot normalized waveforms from the max channel along with template
    normalized = preprocessing.normalize(wf[idx][:,max_chan[idx],:])
    axe[2].plot(normalized.T, lw=0.3, alpha=0.3)
    normalized = normalized_templates[idx]
    axe[2].plot(normalized[max_chan[idx]].T, lw=0.8, color='black')
    axe[0].plot(templates[idx][max_chan[idx]:max_chan[idx] + 1, :].T, lw=0.8)
    # plot normalized template in comparison with other templates
    for i in range(len(wf)):
        if idx == i:
            axe[1].plot(normalized_templates[i][max_chan[i]].T, label=i+1)
        else:
            axe[1].plot(normalized_templates[i][max_chan[i]].T, label=i+1, alpha=0.3)
    # randomly select waveforms from units to plot
    values = randint(0, len(wf), 3)
    for i, val in enumerate(values):    
        axe[3 + i].plot(wf[idx][val, max_chan[idx], :].T)
    # plot pca scores of unit compared with other units and if multiple channel, use max channel for values
    if num_channels > 1:
      for i in range(len(wf)):
          if i != idx:
              axe[6].plot(pca_scores[i][:, max_chan[i], 0], pca_scores[i][:, max_chan[i], 1], 'y*', alpha=0.3)
              axe[7].plot(pca_scores[i][:, max_chan[i], 0], pca_scores[i][:, max_chan[i], 2], 'y*', alpha=0.3)
              axe[8].plot(pca_scores[i][:, max_chan[i], 1], pca_scores[i][:, max_chan[i], 2], 'y*', alpha=0.3)
      axe[6].plot(pca_scores[idx][:, max_chan[idx], 0], pca_scores[idx][:, max_chan[idx], 1], 'r*')
      axe[7].plot(pca_scores[idx][:, max_chan[idx], 0], pca_scores[idx][:, max_chan[idx], 2], 'r*')
      axe[8].plot(pca_scores[idx][:, max_chan[idx], 1], pca_scores[idx][:, max_chan[idx], 2], 'r*')
    else:
      sw.plot_pca_features(recording, sorting, colormap='binary', figure = fig, nproj=3, axes = axe[6:9], max_spikes_per_unit=100)
      sw.plot_pca_features(recording, sorting, unit_ids = [idx+1], figure = fig, nproj=3, axes = axe[6:9] ,max_spikes_per_unit=100)
    # plot isi distribution, amplitude distribution and unit spike train
    sw.plot_isi_distribution(sorting, unit_ids = [idx + 1], bins=100, window=1, axes = axe[9: 11])
    sw.plot_amplitudes_distribution(recording, sorting, max_spikes_per_unit=300, unit_ids = [idx + 1], axes = axe[10: 12])
    sp = sorting.get_unit_spike_train(unit_id = idx + 1)
    axe[11].plot(sp)
    axe[1].legend(loc='lower right', shadow=True, ncol=2)
    axe[0].set_title('Template')
    axe[2].set_title('Waveforms, template')
    axe[3].set_title('Random waveforms')
    axe[6].set_title('PCA component analysis')
    axe[9].set_title('ISI distribution')
    axe[10].set_title('Amplitude distribution')
    axe[11].set_title('Spike frame')
    fig.set_size_inches((8.5, 11), forward=False)
    fig.suptitle("File name: " + file_name + '\nChannel: ' + str(max_chan[idx] + 1).zfill(2) + '\nUnit: ' + str(idx + 1).zfill(2) + '\nFiring Rate: ' + str(round(metrics.loc[idx+1]['firing_rate'], 3)) + ' ISI Violation: ' + str(round(metrics.loc[idx+1]['isi_violation'], 3)) + ' snr: ' + str(round(metrics.loc[idx+1]['snr'], 3)))
    fig.savefig(output_dir + file_name + '_ch' + str(max_chan[idx] + 1).zfill(2) + '_u'+ str(idx+1).zfill(2) + '.png', dpi=500)
    plt.close()
    '''

def maybe_transpose(ts):
  if ts.shape[0] > ts.shape[1]:
    return ts.T
  else:
    return ts

if __name__ == "__main__":
  # directory_in_str = input("Input directory location here")
  directory_in_str = "/media/chang/T3/ml4mat"
  logging.basicConfig(filename='sorting_error.log', level=logging.WARNING)
  directory = os.fsencode(directory_in_str)
  sessions_sorted = []
  if not os.path.exists('pl2_used.txt'):
    open('pl2_used.txt', 'w').close()
  with open('pl2_used.txt', 'r') as f:
    for line in f:
        sessions_sorted.append(line[:-1])
  for file in os.listdir(directory):
    try:
      file_name = os.fsdecode(file)
      file_path = os.path.join(directory_in_str, file_name)
      file_name = file_name.replace('.mat', '')
      print(file_name)
      is_h5_file = True
      if file_name in sessions_sorted:
        continue
      try:
        f = h5py.File(file_path, 'r')
      except:
        f = scipy.io.loadmat(file_path)
        is_h5_file = False  
      ts = np.array(f['mat'])
      ts = maybe_transpose(ts)
      if is_h5_file:
        f.close()
      num_channels = ts.shape[0]
      print('Extracting recording with {} channels ..'.format(num_channels))
      recording = extract_recording(ts, num_channels)
      print('Extracted recording. Preprocessing ..')
      recording_f = recording
      # recording_f = preprocess_recording(recording_f, 300, 6000)
      print('Preprocessed. Sorting ..')
      sorting = sort_recording(recording_f, file_name)
      print('Sorted.')
      wf, max_chan, templates, metrics = postprocess_recording(recording, sorting, file_name)
      output_dir = "waveform_visualization/"
      if not os.path.exists(output_dir):
        os.mkdir(output_dir)
      #visualize_units(recording, recording_f, sorting, wf, templates, max_chan, num_channels, file_name)
      metrics.to_pickle("waveform_visualization/" + file_name + "/metrics_results.pkl")
      with open('pl2_used.txt', 'a') as f:
        f.write('%s\n' % file_name)
      del wf, max_chan, templates
    except Exception as e:
      print('An error occurred processing {}'.format(e))
      logging.exception(file_name)    
    

