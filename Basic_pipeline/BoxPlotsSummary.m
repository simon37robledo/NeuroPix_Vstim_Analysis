%% Do the sumary graphs. Boxplot for moving ball. 

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

examplesSDG =[1 2 3 4 5 6 7 29 30 31 32 40 41 42 43];

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

    goodNeurons = load(sprintf('pvalTime-%s',NP.recordingName)).pvalTi;

    goodNeurons = find(goodNeurons<0.05);

    OSIi = load(sprintf('Orientation-Tuning-Index-%s',NP.recordingName)).L';
    DSIi = load(sprintf('Direction-Selectivity-Index-%s',NP.recordingName)).DSI';

    %theta{i} = load(sprintf('Angle-prefer-%s',NP.recordingName)).Theta;
    OSI{i} = OSIi(goodNeurons);
    DSI{i} = DSIi(goodNeurons);

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
rng(42)
DSIc = [DSI{:};ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))];

mean(ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}])))

OSIc = [OSI{:};ones(1,length([DSI{:}]))+1+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))];

indTOP_DSI= find(DSIc(1,:)>prctile(DSIc(1,:),90));

figure;tiledlayout(1,2,"TileSpacing","compact")

nexttile

% Scatter plot for the first group of points
scatter(DSIc(2,indTOP_DSI), DSIc(1,indTOP_DSI), 50, 'b', 'filled');
hold on;
% Scatter plot for the second group of points
scatter(OSIc(2,indTOP_DSI), OSIc(1,indTOP_DSI), 50, 'r', 'filled');
xticks([mean(ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))),mean(ones(1,length([DSI{:}]))+1+0.1+(0.2-0.1)*rand(1,length([DSI{:}])))])
xticklabels([])

plot([DSIc(2,indTOP_DSI);OSIc(2,indTOP_DSI)],[DSIc(1,indTOP_DSI);OSIc(1,indTOP_DSI)],'Color','k');

ylabel('Score (0-1)');
title('Top 90% percentile DSI');
xlim([1 2.3])
grid on;

% Plot top OSI with lines
rng(42)
DSIc = [DSI{:};ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))];

mean(ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}])))

OSIc = [OSI{:};ones(1,length([DSI{:}]))+1+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))];

indTOP_OSI= find(OSIc(1,:)>prctile(OSIc(1,:),90));

nexttile;

% Scatter plot for the first group of points
scatter(DSIc(2,indTOP_OSI), DSIc(1,indTOP_OSI), 50, 'b', 'filled');
hold on;
% Scatter plot for the second group of points
scatter(OSIc(2,indTOP_OSI), OSIc(1,indTOP_OSI), 50, 'r', 'filled');
xticks([mean(ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))),mean(ones(1,length([DSI{:}]))+1+0.1+(0.2-0.1)*rand(1,length([DSI{:}])))])
xticklabels([])

plot([DSIc(2,indTOP_OSI);OSIc(2,indTOP_OSI)],[DSIc(1,indTOP_OSI);OSIc(1,indTOP_OSI)],'Color','k');
title('Top 90% percentile OSI');
legend('DSI','OSI')
xlim([1 2.3])
grid on;
set(gcf,'Color','w')


%%
% Set the legend
legend('Connecting Line', 'Group 1', 'Group 2', 'Location', 'northwest');

figure;
s = swarmchart(ones(1,length(DSIc(1,indTOP_DSI))),DSIc(1,indTOP_DSI),"filled");
hold on
swarmchart(ones(1,length(DSIc(1,indTOP_DSI)))+1,OSIc(1,indTOP_DSI),"filled");

s.ZJitter = "rand"
s.XJitterWidth = 0.5;

%%
x = [ones(1,500) 2*ones(1,500) 3*ones(1,500)];
y1 = 2 * randn(1,500);
y2 = 3 * randn(1,500) + 5;
y3 = 5 * randn(1,500) + 5;
y = [y1 y2 y3];
s = swarmchart(x,y);




