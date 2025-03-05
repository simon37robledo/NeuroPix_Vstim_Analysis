%%Call NeuronPlot

ex = 51;
cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')

saveDir = '\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure';
eNeuron = 1;



cd(NP.recordingDir)

pvalTi= load(sprintf('pvalsBaselineBoot-1000-%s',NP.recordingName)).pvalsResponse;

respU = find(pvalTi <0.05);

eNeuron = 89;

%% Tuningplot
NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'tuningPlot',1)

%% Raster
NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'raster',1) 

%% Eye positions

NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'EyeMovements',1) 

%% Receptive fields
NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'Refield',1) 
