%%%RF analysis

ex =5;



cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

for ex =  [49 50 52 53 54]
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

    N_bootstrap =1000;

    cd(NP.recordingDir)
    respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
    sign = 0.05; %%%Significance level used to calculate receptive fields
    respU = find(respNeuronsMB<sign);

    NeuronPlotMovingBall(data,ex,respU,...
    'savePlot',1,'saveDir','W:\Large_scale_mapping_NP\Figs paper\receptiveFieldExamples','ReField',1,'noEyeMoves',1,'DivisionType', 'XY')

end
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