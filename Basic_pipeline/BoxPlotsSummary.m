%% Do the sumary graphs. Boxplot for moving ball. 

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

examplesSDG =[28 29 30 31 40 41 42 43];

%%
i=1;
OSI =cell(1,length(examplesSDG));
DSI = cell(1,length(examplesSDG));
for ex = examplesSDG

     path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    try %%In case it is not run in Vstim computer, which has drives mapped differently
        cd(path)
    catch
        originP = cell2mat(extractBetween(path,"\\","\Large_scale"));
        if strcmp(originP,'sil3\data')
            path = replaceBetween(path,"","\Large_scale","W:");
        else
            path = replaceBetween(path,"","\Large_scale","Y:");
        end
        cd(path)
    end
    NP = NPAPRecording(path);

    NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;

    Z_scores = squeeze(NeuronVals(:,:,4));

    rN = find(Z_scores>=prctile(Z_scores,75));


    %theta{i} = load(sprintf('Angle-prefer-%s',NP.recordingName)).Theta;
    OSI{i} = load(sprintf('Orientation-Tuning-Index-%s',NP.recordingName)).L';
    DSI{i} = load(sprintf('Direction-Selectivity-Index-%s',NP.recordingName)).DSI';

    i=i+1;

end
%%

group_names = {'OSI', 'DSI'};

h = figure;
daviolinplot({[OSI{:}]' [DSI{:}]'} ,'boxcolors','k','outliers',0,'whiskers',1,...
    'box',2,'boxwidth',1.2,'scatter',1,'scattersize',10,'jitter',1,'violinwidth',1,'jitterspacing',0.1,...
    'xtlabels', {group_names{1} group_names{2}},'linkline',1);
%title(sprintf('Diff-SDG-%s',strrep(NP.recordingName,'_','-')))
ylabel('Orientation Index - Selectivity Index')
set(h,'Color','w');%yline(-10:2.5:15,'LineWidth',0.1,'Alpha',0.3);%ylim([-10,15]);
cd(NP.recordingDir+"\Figs")
print(h, sprintf('SDG-%s.png',NP.recordingName),'-dpng');

%% Plot top DSI with lines
DSIc = [DSI{:}];

indTOP_DSI= find(DSIc>prctile(DSIc,90));

DSIc(indTOP_DSI);



