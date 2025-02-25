%%%RF analysis

ex =51;

cd(NP.recordingDir)
respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
sign = 0.05;
respU = find(respNeuronsMB<sign);

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);


NeuronPlotMovingBall(data,ex,respU,...
    'savePlot',1,'saveDir','W:\Large_scale_mapping_NP\Figs paper\receptiveFieldExamples','ReField',1,'noEyeMoves',1,'DivisionType', 'XY')


%% %%Plot several examples PV97 & PV35
cd('D:\Mark_S13\Desktop\receptiveFieldExamples');
file = dir ('D:\Mark_S13\Desktop\receptiveFieldExamples');
filenames = {file.name};

append_pdfs(sprintf('RF-%s',NP.recordingName),filenames{:})



% crossCor 
% 
% cc = xcorr2(offsetTemplate,template);
% [max_cc, imax] = max(abs(cc(:)));
% [ypeak, xpeak] = ind2sub(size(cc),imax(1));
% corr_offset = [(ypeak-size(template,1)) (xpeak-size(template,2))];