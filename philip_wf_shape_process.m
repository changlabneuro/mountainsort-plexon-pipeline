clear all
clc
load('meanWFs.mat', 'meanWFs');
% ----------------------------------------------------------------------------------------------------------------------------------------------------
%   Processing of mean WF shape         Processing of mean WF shape         Processing of mean WF shape         Processing of mean WF shape         
%-----------------------------------------------------------------------------------------------------------------------------------------------------
wf_length = 32;
wf_uS = [1:25:800]-201;%[-200:25:575];
interp_wf_length = wf_length*10;
interp_wf_uS = [wf_uS(1): ((wf_uS(end)-wf_uS(1))/interp_wf_length) : wf_uS(end)-((wf_uS(end)-wf_uS(1))/interp_wf_length)];
numWFs = size(meanWFs,1);
wf_matrix = nan(numWFs, wf_length);
norm_wf_matrix = nan(numWFs, wf_length);
interp_wf_matrix = nan(numWFs, interp_wf_length);
trough_peak = nan(numWFs, 1);
half_width = nan(numWFs, 1);
full_amplitude = nan(numWFs, 1);
extremum_amplitude = nan(numWFs, 1);
repolarization_amplitude = nan(numWFs, 1);
for w = 1:numWFs
   
    % Get the analog trace for this unit's mean WF
    this_wf = meanWFs{w};
    
    % Create a normalized amplitude mean WF
    norm_wf = this_wf - min(this_wf(:));
    norm_wf = norm_wf ./ max(norm_wf(:)); % *
    
    % If this waveform was inverted, un-invert it
    if norm_wf(9) == 1
        norm_wf = (-norm_wf)+1;
        this_wf = (-this_wf);
    end
    
    % Create an upsampled interpolated mean WF
    interp_wf = spline(wf_uS,this_wf,interp_wf_uS);
    
    % Save WF's in the matrices
    wf_matrix(w, :)         = this_wf;
    norm_wf_matrix(w, :)    = norm_wf;
    interp_wf_matrix(w, :)  = interp_wf;
    
    above_zero_idxs = find(interp_wf>0);
    below_zero_idxs = find(interp_wf<0);
    
    % Find the post-trough indexs of the interp WF
    interp_wf_postTroughIdxs = find(interp_wf_uS>=0);
    
    % Find the peak and trough of the mean WF
    [~,this_wf_trough_idx] = min(interp_wf);
    [~,this_wf_peak_idx] = max(interp_wf(interp_wf_postTroughIdxs));
    this_wf_peak_idx = this_wf_peak_idx+(interp_wf_postTroughIdxs(1)-1);
    
    % Find the extent of the extremum area
    extremum_start_idx = this_wf_trough_idx;
    extremum_stop_idx = this_wf_trough_idx;
    
    while  extremum_start_idx > 1 && interp_wf(extremum_start_idx-1) <= 0 
        extremum_start_idx = extremum_start_idx-1;
    end
    
    while extremum_stop_idx <= interp_wf_length-1 && interp_wf(extremum_stop_idx-1) <= 0 
        extremum_stop_idx = extremum_stop_idx+1;
    end
    
    extremum_X = interp_wf_uS(extremum_start_idx:extremum_stop_idx);
    extremum_Y = interp_wf(extremum_start_idx:extremum_stop_idx);
    extremum_Q = trapz(extremum_X,extremum_Y);
    
     % Find the extent of the repolarization area
    repolarization_start_idx = this_wf_peak_idx;
    repolarization_stop_idx = this_wf_peak_idx;
     while  repolarization_start_idx > 1 && interp_wf(repolarization_start_idx-1) >= 0 
        repolarization_start_idx = repolarization_start_idx-1;
    end
    
    while repolarization_stop_idx <= interp_wf_length-1 && interp_wf(repolarization_stop_idx-1) >= 0 
        repolarization_stop_idx = repolarization_stop_idx+1;
    end
     
    repolarization_X = interp_wf_uS(repolarization_start_idx:repolarization_stop_idx);
    repolarization_Y = interp_wf(repolarization_start_idx:repolarization_stop_idx);
    repolarization_Q = trapz(repolarization_X,repolarization_Y);
      
      
     % Find the extent of the capacitance
     
     if extremum_start_idx > 1
         
         capacitance_stop_idx = extremum_start_idx;
         capacitance_start_idx = extremum_start_idx;
         
           while capacitance_start_idx > 1 && interp_wf(capacitance_start_idx-1) >= 0 
               capacitance_start_idx = capacitance_start_idx-1;
           end
           capacitance_X = interp_wf_uS(capacitance_start_idx:capacitance_stop_idx);
            capacitance_Y = interp_wf(capacitance_start_idx:capacitance_stop_idx);
            capacitance_Q = trapz(capacitance_X,capacitance_Y);
         
     else
         capacitance_X = [nan];
         capacitance_Y = [nan];
         capacitance_Q = 0;
     end
    
        extremum_auc = extremum_Q;
    repolarization_auc = repolarization_Q;
    
    
    this_wf_peak_uS = interp_wf_uS(this_wf_peak_idx);
    this_wf_trough_uS = interp_wf_uS(this_wf_trough_idx);
    
    fprintf('Peak: %10.5f�S\n', this_wf_peak_uS);
    fprintf('Trough: %10.5f�S\n', this_wf_trough_uS);
    
    this_wf_t2p_uS = this_wf_peak_uS-this_wf_trough_uS;
    
    fprintf('Trough-to-peak: %10.5f�S\n', this_wf_t2p_uS);
   
    wf_trough2peak_uS(w) = this_wf_t2p_uS;
    
    
    
    extremum_amp = interp_wf(this_wf_trough_idx);
    repolarization_amp = interp_wf(this_wf_peak_idx);
    
    fprintf('Extremum amplitude: %7.5f�V\t/\tAUC: %8.5f\n', extremum_amp, extremum_Q);
    fprintf('Repolarization amplitude: %7.5f�V\t/\tAUC: %8.5f\n', repolarization_amp, repolarization_Q);
    
    half_width_point_amp = extremum_amp/2;
    
    
    decending_extremum_idxs = extremum_start_idx:this_wf_trough_idx;
    decending_extremum_X = interp_wf_uS(decending_extremum_idxs);
    decending_extremum_Y = interp_wf(decending_extremum_idxs);
    
    [~, decending_extremum_half_width_idx] = min(abs(decending_extremum_Y-half_width_point_amp));
    decending_extremum_half_width_X = decending_extremum_X(decending_extremum_half_width_idx);
    decending_extremum_half_width_Y = decending_extremum_Y(decending_extremum_half_width_idx);
    
    accending_extremum_idxs =  this_wf_trough_idx:extremum_stop_idx;
    accending_extremum_X = interp_wf_uS(accending_extremum_idxs);
    accending_extremum_Y = interp_wf(accending_extremum_idxs);
    
    [~, accending_extremum_half_width_idx] = min(abs(accending_extremum_Y-half_width_point_amp));
    accending_extremum_half_width_X = accending_extremum_X(accending_extremum_half_width_idx);
    accending_extremum_half_width_Y = accending_extremum_Y(accending_extremum_half_width_idx);
    
    half_width_point_uS = accending_extremum_half_width_X-decending_extremum_half_width_X;
    
    
    figure(1);
    clf;
    hold on;
    
    plot(wf_uS, this_wf, 'LineWidth', 2);
    plot(interp_wf_uS, interp_wf, 'LineWidth', 2);
    hline(0, 'k--');
    
  
    
    plot([this_wf_peak_uS this_wf_peak_uS], [0 repolarization_amp], 'g--', 'LineWidth', 1.5)
    
    plot([this_wf_trough_uS this_wf_trough_uS], [0 extremum_amp], 'c--', 'LineWidth', 1.5)
    
      patch('XData',extremum_X,'YData',extremum_Y,'FaceColor','cyan','EdgeColor','none', 'FaceAlpha', 0.1)
      patch('XData',repolarization_X,'YData',repolarization_Y,'FaceColor','green','EdgeColor','none', 'FaceAlpha', 0.1)
      
       patch('XData',capacitance_X,'YData',capacitance_Y,'FaceColor','yellow','EdgeColor','none', 'FaceAlpha', 0.1)
      
       plot([decending_extremum_half_width_X accending_extremum_half_width_X], [decending_extremum_half_width_Y accending_extremum_half_width_Y], 'LineWidth', 2);
      
    pause
    
    trough_peak(w) = this_wf_t2p_uS;
    half_width(w) = half_width_point_uS;
    full_amplitude(w) = abs(extremum_amp)+abs(repolarization_amp);
    extremum_amplitude(w) = extremum_amp;
    repolarization_amplitude(w) = repolarization_amp;
    
    clc
end
clust_1_color = [0                     0.447                     0.741];
clust_2_color = [0.85                     0.325                     0.098];
pca_matrix = [trough_peak, half_width, full_amplitude, extremum_amplitude, repolarization_amplitude];
[pca_coeff,pca_score,pca_latent] = pca(pca_matrix);
X = interp_wf_matrix;
Y = pca_score;
clust_idxs = kmeans(Y,2);
clust_1_idxs = find(clust_idxs == 1);
clust_2_idxs = find(clust_idxs == 2);
plotH = figure(1);
subplot(3,2,3);
hold on;
hist_edges = [min(Y(:,1)) : ((max(Y(:,1))-min(Y(:,1)))/24) : max(Y(:,1))];
histogram(Y(clust_1_idxs,1), hist_edges, 'EdgeColor', 'none', 'FaceColor',clust_1_color)
histogram(Y(clust_2_idxs,1), hist_edges, 'EdgeColor', 'none', 'FaceColor',clust_2_color)
xlabel('1st component')
ylabel('Count (# of units)')
title('Histogram of 1st PCA component');
% subplot(3,2,2);
% hold on;
% hist_edges = [min(Y(:,2)) : ((max(Y(:,2))-min(Y(:,2)))/24) : max(Y(:,2))];
% histogram(Y(clust_1_idxs,2), hist_edges, 'EdgeColor', 'none', 'FaceColor',clust_1_color)
% histogram(Y(clust_2_idxs,2), hist_edges, 'EdgeColor', 'none', 'FaceColor',clust_2_color)
% xlabel('2nd component')
% ylabel('Count (# of units)')
% title('Histogram of 2nd PCA component');
subplot(3,2,1);
hold on;
title('PCA with k-means clustering (k=2)');
scatter(Y(clust_1_idxs,1),Y(clust_1_idxs,2),55,clust_1_color,'filled');
scatter(Y(clust_2_idxs,1),Y(clust_2_idxs,2),55,clust_2_color,'filled');
xlabel('1st component')
ylabel('2nd component')
subplot(3,2,4);
hold on;
clust_1_wfs = X(clust_1_idxs, :);
clust_2_wfs = X(clust_2_idxs, :);
clust_1_mean_wf = mean(clust_1_wfs);
clust_2_mean_wf = mean(clust_2_wfs);
plot(interp_wf_uS, clust_1_mean_wf, 'Color', clust_1_color, 'LineWidth', 3);
plot(interp_wf_uS, clust_2_mean_wf, 'Color', clust_2_color, 'LineWidth', 2);
xlabel('\muS')
ylabel('Upsampled WF amplitude (+/- 1 STD)')
clust_1_std_wf = std(clust_1_wfs);
clust_2_std_wf = std(clust_2_wfs);
clust_1_ste_wf = std(clust_1_wfs)/sqrt(size(clust_1_idxs,1));
clust_2_ste_wf = std(clust_2_wfs)/sqrt(size(clust_2_idxs,1));
clust_1_upperVar = [clust_1_mean_wf+(clust_1_std_wf)];
clust_1_lowerVar =[clust_1_mean_wf-(clust_1_std_wf)];
clust_1_patchX = horzcat([interp_wf_uS], fliplr([interp_wf_uS]));
clust_1_patchY = horzcat(clust_1_upperVar, fliplr(clust_1_lowerVar));
clust_1_H = patch('XData', clust_1_patchX,'YData', clust_1_patchY,'FaceColor',clust_1_color, 'FaceAlpha',0.3, 'EdgeColor', 'none');
clust_2_upperVar = [clust_2_mean_wf+(clust_2_std_wf)];
clust_2_lowerVar =[clust_2_mean_wf-(clust_2_std_wf)];
clust_2_patchX = horzcat([interp_wf_uS], fliplr([interp_wf_uS]));
clust_2_patchY = horzcat(clust_2_upperVar, fliplr(clust_2_lowerVar));
clust_2_H = patch('XData', clust_2_patchX,'YData', clust_2_patchY,'FaceColor',clust_2_color, 'FaceAlpha',0.3, 'EdgeColor', 'none');
subplot(3,2,5);
hold on;
title('Half-width');
clust_1_half_widths = half_width(clust_1_idxs);
clust_2_half_widths = half_width(clust_2_idxs);
ylabel('\muS')
xlabel('Cluster');
xticks([1 2]);
bplot(clust_1_half_widths,1, 'color', clust_1_color);
bplot(clust_2_half_widths,2, 'color', clust_2_color);
subplot(3,2,6);
hold on;
title('Trough-to-peak');
clust_1_t2p = trough_peak(clust_1_idxs);
clust_2_t2p = trough_peak(clust_2_idxs);
ylabel('\muS')
xlabel('Cluster');
xticks([1 2]);
bplot(clust_1_t2p,1, 'color', clust_1_color);
bplot(clust_2_t2p,2, 'color', clust_2_color);
subplot(3,2,2);
hold on;
title('Amplitude');
clust_1_amp = full_amplitude(clust_1_idxs);
clust_2_amp = full_amplitude(clust_2_idxs);
ylabel('Voltage')
xlabel('Cluster');
xticks([1 2]);
bplot(clust_1_amp,1, 'color', clust_1_color);
bplot(clust_2_amp,2, 'color', clust_2_color);
% 
%  saveStr = sprintf('WF_shape.pdf');
%         savePath = fullfile(pwd, saveStr);
%           set(plotH, 'position', [-1880         10        1800         1000]);
%         
%         set(plotH,'PaperOrientation','landscape');
%         set(plotH,'PaperUnits','normalized');
%         set(plotH,'PaperPosition', [0 0 1 1]);
%         print(plotH, '-dpdf', savePath);