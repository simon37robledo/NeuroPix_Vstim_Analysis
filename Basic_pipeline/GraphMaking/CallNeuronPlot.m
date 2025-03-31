%%Call NeuronPlot

ex = 51;

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

saveDir = '\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure';
cd(saveDir)

NP = loadNPclassFromTable(ex);

cd(NP.recordingDir)

pvals= load(sprintf('pvalsBaselineBoot-1000-%s',NP.recordingName)).pvalsResponse;
sigma = 0.05;
goodNeurons =  find(pvals <sigma);
p = NP.convertPhySorting2tIc(NP.recordingDir);
label = string(p.label');
goodU = p.ic(:,label == 'good');
ampsAll = p.neuronAmp(label == 'good');
phyID = p.phy_ID(label == 'good');

ampsGood = ampsAll(goodNeurons);
find(ampsGood>50)

eNeuron = goodNeurons(82); %%%Neuron for paper

ampsGood(82)

phyID(eNeuron)

ch = goodU(1,eNeuron);

possibleUnits = sort(phyID(goodNeurons(ampsGood>50)));

eNeuron = find(phyID == 10);


%% Tuningplot
NeuronPlotMovingBall(data,ex,eNeuron,...
    'tuningPlot',1,'savePlot',1,'saveDir',saveDir)

%% Raster
NeuronPlotMovingBall(data,ex,eNeuron,...
    'raster',1,'savePlot',1,'saveDir',saveDir) 

%% Eye positions

NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'EyeMovements',1) 

cd(saveDir)
print(gcf, 'EyeMovsPV97_3', '-dpdf', '-r300', '-vector');
     

%% Receptive fields no eye movements
NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'Refield',1,'noEyeMoves',1) 

%% Receptive fields no eye movements
NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'Refield',1) 
