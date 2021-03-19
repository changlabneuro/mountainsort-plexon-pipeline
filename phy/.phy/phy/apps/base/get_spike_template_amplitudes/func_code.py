# first line: 490
    def get_spike_template_amplitudes(self, spike_ids, **kwargs):
        """Return the spike template amplitudes as stored in `amplitudes.npy`."""
        if self.model.amplitudes is None:
            return np.zeros(len(spike_ids))
        amplitudes = self.model.amplitudes[spike_ids]
        return amplitudes
