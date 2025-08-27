cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

ex =82; %error with 67 PV132_4
NP = loadNPclassFromTable(ex); %73 81

vsR = rectNoiseGridAnalysis(NP);
vsR.getSessionTime("overwrite",true);
d2 = vsR.getDiodeTriggers(); 
d3 = vsR.getSyncedDiodeTriggers();

r = vsR.getReceptiveFields();
vsR.plotReceptiveFields
%% Rect Grid
NP = loadNPclassFromTable(84); %73 81
vsRe = rectGridAnalysis(NP);
vsRe.getSessionTime("overwrite",true);
%vsRe.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true); 
vsRe.getDiodeTriggers('overwrite',true); 
vsRe.getSyncedDiodeTriggers("overwrite",true);
vsRe.plotSpatialTuningSpikes;
vsRe.plotSpatialTuningLFP;
%% Moving ball
NP = loadNPclassFromTable(86); %73 81
vs = linearlyMovingBallAnalysis(NP);
vs.getSessionTime("overwrite",true);
vs.getDiodeTriggers('extractionMethod','digitalTriggerDiode','overwrite',true); 
%vs.plotDiodeTriggers
vs.getSyncedDiodeTriggers("overwrite",true);
%vs.plotSpatialTuningSpikes;
r = vs.ResponseWindow('overwrite',true);
results = vs.ShufflingAnalysis('overwrite',true);
vs.plotRaster('AllResponsiveNeurons',true,'overwrite',true,'MergeNtrials',3)
vs.plotRaster('exNeurons',42,'overwrite',true,'MergeNtrials',3)
vs.plotCorrSpikePattern
%vs.plotRaster('AllSomaticNeurons',true,'overwrite',true,'speed',2)
vs.CalculateReceptiveFields('overwrite',true);
vs.PlotReceptiveFields("allDir",true,'overwrite',true)

%% Check experiments in timseseries viewer
timeSeriesViewer(NP)
t=NP.getTrigger;
data.VS_ordered(ex)

stimOn = t{3};
stimOff = t{4};

MBRtOn = stimOn(stimOn > t{1}(1) & stimOn < t{2}(1));
MBRtOff = stimOff(stimOff > t{1}(1) & stimOff < t{2}(1));

MBtOn = stimOn(stimOn > t{1}(2) & stimOn < t{2}(2));
MBtOff = stimOff(stimOff > t{1}(2) & stimOff < t{2}(2));

RGtOn = stimOn(stimOn > t{1}(3) & stimOn < t{2}(3));
RGtOff = stimOff(stimOff > t{1}(3) & stimOff < t{2}(3));

NGtOn = stimOn(stimOn > t{1}(4) & stimOn < t{2}(4));
NGtOff = stimOff(stimOff > t{1}(4) & stimOff < t{2}(4));

DtOn = stimOn(stimOn > t{1}(5) & stimOn < t{2}(5));
DtOff = stimOff(stimOff > t{1}(5) & stimOff < t{2}(5));

MovingBallTriggersDiode = d3.stimOnFlipTimes;



%% %% check neural data sync and analog data sync 

allTimes = [stimOn(:); stimOff(:); onSync(:); offSync(:)];  % concatenate as column

% Sort from earliest to latest
sortedTimesDiodeOldMethod = sort(allTimes);
