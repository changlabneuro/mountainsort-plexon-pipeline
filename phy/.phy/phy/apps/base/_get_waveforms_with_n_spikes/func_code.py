# first line: 182
    def _get_waveforms_with_n_spikes(
            self, cluster_id, n_spikes_waveforms, batch_size_waveforms, current_filter=None):
        # HACK: we pass self.raw_data_filter.current_filter so that it is cached properly.
        pos = self.model.channel_positions
        spike_ids = self.selector.select_spikes(
            [cluster_id], n_spikes_waveforms, batch_size_waveforms)
        channel_ids = self.get_best_channels(cluster_id)
        channel_labels = self._get_channel_labels(channel_ids)
        # # Get all spikes from the cluster, will be selected below depending on the chunks.
        # spike_ids = self.supervisor.clustering.spikes_per_cluster[cluster_id]
        # spike_times = self.model.spike_times[spike_ids]
        # # Subselection to minimize the number of chunks to load.
        # spike_ids = spike_ids[select_spikes_from_chunked(
        #     spike_times, self.model.traces.chunk_bounds, n_spikes_waveforms)]
        data = self.model.get_waveforms(spike_ids, channel_ids)
        if isinstance(data, da.Array):  # pragma: no cover
            data = data.compute()
        assert data.ndim == 3  # n_spikes, n_samples, n_channels
        # data = data - data.mean() if data is not None else None
        # Filter the waveforms.
        if data is not None:
            data = self.raw_data_filter.apply(data, axis=1)
        return Bunch(
            data=data,
            channel_ids=channel_ids,
            channel_labels=channel_labels,
            channel_positions=pos[channel_ids],
        )
