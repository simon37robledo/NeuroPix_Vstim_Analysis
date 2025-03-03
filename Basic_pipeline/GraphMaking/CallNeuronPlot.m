%%Call NeuronPlot

ex = 51;
saveDir = 'W:\Large_scale_mapping_NP\Figs paper\1stFigure';
eNeuron = 14;

%% Tuningplot
NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'tuningPlot',1)

%% Raster
NeuronPlotMovingBall(data,ex,eNeuron,...
    'savePlot',1,'saveDir',saveDir,'raster',1) 
