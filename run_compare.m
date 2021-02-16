old_spiketimes = load('dictator_game_SUAdata_pre.mat');
new_spiketimes = load('spiketimes.mat');
compare_out = util.compare_spiketimes(old_spiketimes, new_spiketimes);