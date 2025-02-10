%% ~~ Example bombcell pipeline ~~
% Adjust the paths in the 'set paths' section and the parameters in bc_qualityParamValues
% This pipeline will:
%   (1) load your ephys data, 
%   (2) decompress your raw data if it is in .cbin format 
%   (3) run bombcell on your data and save the output and
%   (4) bring up summary plots and a GUI to flip through classified cells.
% The first time, this pipeline will be significantly slower (10-20' more)
% than after because it extracts raw waveforms. Subsequent times these
% pre-extracted waveforms are simply loaded in.
% We recommend running this pipeline on a few datasets and deciding on
% quality metric thresholds depending on the summary plots (histograms 
% of the distributions of quality metrics for each unit) and GUI. 


%% set paths - EDIT THESE 
% ephysKilosortPath = '\\132.66.45.127\data\\Large_scale_mapping_NP\Awake_exp\SA6\SA6_experiment_23_8_23\Insertion1\catgt_SA6_experiment_23_8_23_1_g0\';% path to your kilosort output files 
% ephysRawDir = dir('\\132.66.45.127\data\\Large_scale_mapping_NP\Awake_exp\SA6\SA6_experiment_23_8_23\Insertion1\catgt_SA6_experiment_23_8_23_1_g0\SA6_experiment_23_8_23_1_g0_tcat.imec0.ap.bin'); % path to yourraw .bin or .dat data
% ephysMetaDir = dir('\\132.66.45.127\data\\Large_scale_mapping_NP\Awake_exp\SA6\SA6_experiment_23_8_23\Insertion1\catgt_SA6_experiment_23_8_23_1_g0\SA6_experiment_23_8_23_1_g0_tcat.imec0.ap.meta'); % path to your .meta or .oebin meta file
% saveLocation = '\\132.66.45.127\data\Large_scale_mapping_NP\Awake_exp\SA6\SA6_experiment_23_8_23\Insertion1\QM_bombcell'; % where you want to save the quality metrics 
% savePath = fullfile(saveLocation, 'qMetrics'); 
% decompressDataLocal = '\\132.66.45.127\data\Large_scale_mapping_NP\Awake_exp\SA6\SA6_experiment_23_8_23\Insertion1\QM_bombcell'; % where to save raw decompressed ephys data 

ephysKilosortPath = '\\132.66.45.127\data\\Large_scale_mapping_NP\Immobilized_exp\PV67\PV67_experiment_5_7_23\Insertion1\catgt_PV67_Experiment_5_7_23_1_g0\';% path to your kilosort output files 
ephysRawDir = dir('\\132.66.45.127\data\\Large_scale_mapping_NP\Immobilized_exp\PV67\PV67_experiment_5_7_23\Insertion1\catgt_PV67_Experiment_5_7_23_1_g0\PV67_Experiment_5_7_23_1_g0_tcat.imec0.ap.bin'); % path to yourraw .bin or .dat data
ephysMetaDir = dir('\\132.66.45.127\data\\Large_scale_mapping_NP\Immobilized_exp\PV67\PV67_experiment_5_7_23\Insertion1\catgt_PV67_Experiment_5_7_23_1_g0\PV67_Experiment_5_7_23_1_g0_tcat.imec0.ap.meta'); % path to your .meta or .oebin meta file
saveLocation = '\\132.66.45.127\data\Large_scale_mapping_NP\Immobilized_exp\PV67\PV67_experiment_5_7_23\Insertion1'; % where you want to save the quality metrics 
savePath = fullfile(saveLocation, 'qMetrics'); 
decompressDataLocal = '\\132.66.45.127\data\Large_scale_mapping_NP\Awake_exp\SA6\SA6_experiment_23_8_23\Insertion1'; % where to save raw decompressed ephys data 

%% load data 
[spikeTimes_samples, spikeTemplates, templateWaveforms, templateAmplitudes, pcFeatures, ...
    pcFeatureIdx, channelPositions] = bc_loadEphysData(ephysKilosortPath);

%% detect whether data is compressed, decompress locally if necessary
rawFile = bc_manageDataCompression(ephysRawDir, decompressDataLocal);

%% which quality metric parameters to extract and thresholds 
param = bc_qualityParamValues(ephysMetaDir, rawFile, ephysKilosortPath); %for unitmatch, run this:
% param = bc_qualityParamValuesForUnitMatch(ephysMetaDir, rawFile, ephysKilosortPath, gain_to_uV)

%% compute quality metrics 
rerun = 0;
qMetricsExist = ~isempty(dir(fullfile(savePath, 'qMetric*.mat'))) || ~isempty(dir(fullfile(savePath, 'templates._bc_qMetrics.parquet')));

if qMetricsExist == 0 || rerun
    [qMetric, unitType] = bc_runAllQualityMetrics(param, spikeTimes_samples, spikeTemplates, ...
        templateWaveforms, templateAmplitudes, pcFeatures, pcFeatureIdx, channelPositions, savePath);
else
    [param, qMetric] = bc_loadSavedMetrics(savePath); 
    unitType = bc_getQualityUnitType(param, sortrows(qMetric,'maxChannels'), savePath);
end


%% view units + quality metrics in GUI 
% load data for GUI
loadRawTraces = 1; % default: don't load in raw data (this makes the GUI significantly faster)
bc_loadMetricsForGUI;

% GUI guide: 
% left/right arrow: toggle between units 
% g : go to next good unit 
% m : go to next multi-unit 
% n : go to next noise unit
% up/down arrow: toggle between time chunks in the raw data
% u: brings up a input dialog to enter the unit you want to go to
unitQualityGuiHandle = bc_unitQualityGUI(memMapData, ephysData, qMetric, forGUI, rawWaveforms, ...
    param, probeLocation, unitType, loadRawTraces);
%% Save figures

cd('\\132.66.45.127\data\\Large_scale_mapping_NP\Immobilized_exp\PV67\PV67_experiment_5_7_23\Insertion1\qMetrics')

savefig(figure(1), 'templateWaveform.fig')
savefig(figure(2), 'qualityMetrics.fig')

close all
%% Asigning spikes to figures

%0; % Noise
%1; % Good Somatic
%2; % MUA
%3; % NON-Somatic
cd('\\132.66.45.127\data\\Large_scale_mapping_NP\Immobilized_exp\PV67\PV67_experiment_5_7_23\Insertion1\qMetrics')
path = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV67\PV67_experiment_5_7_23\Insertion1\catgt_PV67_experiment_5_7_23_1_g0';

NP = NPAPRecording(path);

qMetric=sortrows(qMetric,'maxChannels');

dspikes = readNPY('spikes._bc_duplicateSpikes.npy');

GoodUtemplate = find(unitType ~= 0);

ic = zeros(5,numel(GoodUtemplate)); %cell(2, length(GoodUtemplate));

nonDuplicateSPKs = find(dspikes == 0);

spikeTemplatesND = spikeTemplates(nonDuplicateSPKs);

spikeTimes_samplesND = spikeTimes_samples(nonDuplicateSPKs);

t = cell(1,numel(GoodUtemplate));
currentIdx=0;prevCh=-1;
savei=[];

for i = 1:numel(GoodUtemplate)

    if ~isempty(spikeTimes_samplesND(spikeTemplatesND == GoodUtemplate(i)))

        tm = spikeTimes_samplesND(spikeTemplatesND == GoodUtemplate(i));

        t{i} = uint64(tm');

        ic(1,i) = qMetric.maxChannels(GoodUtemplate(i),:);
        ic(3,i)=currentIdx+1;
        ic(4,i)=currentIdx+numel(t{i});
        ic(5,i)=qMetric.phy_clusterID(GoodUtemplate(i),:);

        if prevCh==ic(1,i)
            ic(2,i)=ic(2,i-1)+1;
        else
            ic(2,i)=1;
        end
        prevCh=ic(1,i);
        currentIdx=ic(4,i);

        savei = [savei i];
        
    end
end

t=double(cell2mat(t(find(~cellfun(@isempty,t)))))/(NP.samplingFrequency(1)/1000);
ic = ic(:,ic(1,:)~=0);




%%

 pBC = NP.convertBCSorting2tIc(NP.recordingDir);

 clusterTable=readtable([ephysKilosortPath 'cluster_info.tsv'],'FileType','delimitedtext');

 p = NP.convertPhySorting2tIc(NP.recordingDir);

 %% 
cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);


 %%
allGoodRec = [1:21, 28:36, 40:54];
% %path = '\\sil3\data\\Large_scale_mapping_NP\lizards\PV139\PV139_Experiment_6_2_24\Insertion1\catgt_PV139_Experiment_6_2_24_1_g0';
% path = '\\132.66.45.127\data\Large_scale_mapping_NP\Mice_experiments\mouse1\Mice_exp_28_11_23\Insertion1\catgt_Mice_exp_28_11_23_1_g0';

% Iterate through experiments (insertions and animals) in excel file
for ex =  54%allGoodRec %GoodRecordings%GoodRecordingsPV%GoodRecordingsPV%selecN{1}(1,:) %1:size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
    path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    try %%In case it is not run in Vstim computer, which has drives mapped differently
        cd(pathE)
    catch
        try
            originP = cell2mat(extractBetween(path,"\\","\Large_scale"));
            if strcmp(originP,'sil3\data')
                path = replaceBetween(path,"","\Large_scale","W:");
            else
                path = replaceBetween(path,"","\Large_scale","Y:");
            end

            cd(path)
        catch

            if strcmp(originP,'sil3\data')
                path = replaceBetween(path,"","\Large_scale","\\sil3\data");
            else
                path = replaceBetween(path,"","\Large_scale","\\sil1\data");
            end
            cd(path)

        end

    end

    NP = NPAPRecording(path);

    [qMetric,unitType]=NP.getBombCell(NP.recordingDir,1);
    close all

end




