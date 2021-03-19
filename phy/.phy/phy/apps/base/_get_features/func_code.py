# first line: 357
    def _get_features(self, cluster_id=None, channel_ids=None, load_all=False):
        """Return the features of a given cluster on specified channels."""
        spike_ids = self._get_feature_view_spike_ids(cluster_id, load_all=load_all)
        # Use the best channels only if a cluster is specified and
        # channels are not specified.
        if cluster_id is not None and channel_ids is None:
            channel_ids = self.get_best_channels(cluster_id)
        return self._get_spike_features(spike_ids, channel_ids)
