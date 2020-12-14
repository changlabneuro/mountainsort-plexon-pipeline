function [success, mean_wf] = extractMeanWF(uniqueID, tablePath)
%EXTRACTMEANWF Summary of this function goes here
%   Detailed explanation goes here
% -----------------------	TEST SETTINGS   -----------------------
% clear all
% clc
% tablePath = 'C:\Users\Bread\Desktop\combinedTable\nestedTable2.mat';
% uniqueID = 'Mittens_112916_MCR1_C3U1';%'Ham_101217_FIO11_C30U1';%'Ham_101017_FIO10_C19U2';
% ----------------------- Positive thinking -----------------------
success = true;
extractNewAnalog = false;
% ----------------------- Settings -----------------------
preThresholdUs = 600;
postThresholdUs = 1400;
preThresholdS = preThresholdUs/1e6;
postThresholdS = postThresholdUs/1e6;
% -------------- Load table and locate our desired unit  --------------
loadStruct = load(tablePath);
existingTable = loadStruct.cellTable;
unit_TableIdx = find(strcmpi(existingTable.UniqueID, uniqueID));
cellStruct = table2struct(existingTable);
if isempty(unit_TableIdx)
    error('Desired unit not found in table!!');
end
% -------------- Load the spike times --------------
try
    spikeFilePath = cellStruct(unit_TableIdx).SpikePath;
    load(spikeFilePath, 'spikes');
catch spikeLoadErr
    success = false;
    warning(['XNT (add2Table_StimOn_Response): Could not load spike times or waveforms --> ', spikeLoadErr.message]);
    return
end
% Extract the session and unit strings from the table
unitStr = cellStruct(unit_TableIdx).Unit;
sessionStr = cellStruct(unit_TableIdx).Session;
dataPath = cellStruct(unit_TableIdx).DataPath;
validSpikeLookups = find(~cellfun(@isempty,spikes));
unitNames = cell(size(validSpikeLookups));
% --------------- Loop through valid units to find our desired unit ---------------
for i = 1:size(validSpikeLookups,1)
    
    u = validSpikeLookups(i);
    
    unitStruct = spikes{u};
    
    if isfield(unitStruct, 'trialLimits')
        unitLimits = unitStruct.trialLimits;
    else
        unitLimits = [];
    end
    
    unitName = unitStruct.name;
    
    unitNames{i} = unitName;
    
end
desiredUnit_Idx = find(strcmp(unitNames, unitStr));
desiredUnit_Struct =  spikes{validSpikeLookups(desiredUnit_Idx)};
desiredUnit_Spikes = desiredUnit_Struct.ts;
desiredUnit_SourceChan =  desiredUnit_Struct.infoValues{    find(strcmp(desiredUnit_Struct.infoNames, 'Channel Name')) };
desiredUnit_SourceChan = strrep(desiredUnit_SourceChan, 'SPK_', '');
desiredUnit_SourceChan = strrep(desiredUnit_SourceChan, 'CutFromSMRX', '');
desiredUnit_numSpikes = size(desiredUnit_Spikes,1);
if extractNewAnalog
    
    spikeChannel_Struct = PL2Ad(dataPath, desiredUnit_SourceChan);
    spikeChannel_fs = spikeChannel_Struct.ADFreq;
    spikeChannel_numSamples = spikeChannel_Struct.FragCounts;
    spikeChannel_timeInterval = 1/spikeChannel_fs;
    spikeChannel_maxS = spikeChannel_timeInterval*spikeChannel_numSamples;
    spikeChannel_tv = [0:spikeChannel_timeInterval:(spikeChannel_maxS-spikeChannel_timeInterval)];
    
    extractedAnalog_numSamples = size([(-preThresholdS):spikeChannel_timeInterval:postThresholdS],2);
    
    extractedAnalog = nan(desiredUnit_numSpikes, extractedAnalog_numSamples);
    
    
    for s = 1:desiredUnit_numSpikes
        
        thisSpike_thresholdS = desiredUnit_Spikes(s);
        
        thisSpike_earlierS = thisSpike_thresholdS-preThresholdS;
        thisSpike_laterS = thisSpike_thresholdS+postThresholdS;
        
        if thisSpike_earlierS > 0 && thisSpike_laterS <= spikeChannel_maxS
            
            %This method is apparently slower
            %         tic
            %         [~, thisSpike_earlier_TV_idx] = min(abs(spikeChannel_tv-thisSpike_earlierS));
            %
            %         thisSpike_tvIdxs = [thisSpike_earlier_TV_idx : (thisSpike_earlier_TV_idx+(extractedAnalog_numSamples-1))];
            %         toc
            
            
            
            thisSpike_tvIdxs = find(spikeChannel_tv>= thisSpike_earlierS & spikeChannel_tv <= thisSpike_laterS);
            % Since sometimes we're getting off by one sample let's just take the N samples after the first one for a fixed length vector
            thisSpike_tvIdxs = thisSpike_tvIdxs(1):(thisSpike_tvIdxs(1)+(extractedAnalog_numSamples-1));
            
            
            
            thisSpike_analog =  spikeChannel_Struct.Values(thisSpike_tvIdxs);
            extractedAnalog(s,:) =thisSpike_analog;
            
        end
        
    end
    
    mean_wf = nanmean(extractedAnalog);
elseif isfield(desiredUnit_Struct, 'wf')
    
    all_wfs = desiredUnit_Struct.wf;
    mean_wf = mean(all_wfs);
else
    success = false;
    mean_wf = [];
end