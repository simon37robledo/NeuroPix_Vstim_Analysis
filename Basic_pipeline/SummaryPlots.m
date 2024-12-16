%% Do the sumary graphs. Boxplot for moving ball. 

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

GoodRecordings =[1:20,40:48];%[1:20,28:32,40:48]; Perform bootstrap with SA5

%%
animalN = 7;
i=1;
OSI =cell(1,length(GoodRecordings));
DSI = cell(1,length(GoodRecordings));
PreferAngle = cell(1,length(GoodRecordings));
animalID = cell(1,animalN);
DepthUnit = cell(1,length(GoodRecordings));
RespMB = cell(1,length(GoodRecordings));
RespRG = cell(1,length(GoodRecordings));
goodNeuronsMB =  cell(1,length(GoodRecordings));
NeuronN = cell(1,length(GoodRecordings));
tablePos = cell(1,length(GoodRecordings));

pvalsMB = cell(1,length(GoodRecordings));

pvalsRG = cell(1,length(GoodRecordings));

N_bootstrap = 1000;

for ex = GoodRecordings

     path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    try %%In case it is not run in Vstim computer, which has drives mapped differently
        cd(path)
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


    if contains(data.VS_ordered(ex),"MB")
        NeuronValsMB = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;

        respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        respNeuronsRG = load(sprintf('RectGrid-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        goodNeuronsMB{i} = find(respNeuronsMB<0.005);
        tuningCurve = load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve;
        uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
        [preferDir dirInd] = max(tuningCurve,[],2);
        PreferAngle{i} = preferDir(respNeuronsMB<0.005);
        RespMB{i} = max(NeuronValsMB(respNeuronsMB<0.005,:,4),[],2)';
        NeuronN{i} = find(respNeuronsMB<0.005); 
        %RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
        pvalsMB{i} = respNeuronsMB(respNeuronsMB<0.005);


    else
        goodNeuronsMB{i} =0;
        PreferAngle{i} = 0;
        RespMB{i} = 0;
        NeuronN{i} =0;
        pvalsMB{i} =0;

    end


    if contains(data.VS_ordered(ex),"RG") &&  contains(data.VS_ordered(ex),"MB")

          NeuronValsRG= load(sprintf('RectGrid-NeuronRespCat-%s',NP.recordingName)).NeuronVals;
          %EnoughTrials = NeuronValsRG(respNeuronsMB<0.005,:,[4 7]);

          maxRG = (max(NeuronValsRG(respNeuronsMB<0.005,:,[4 7]),[],2));

          %RespRG{i} = (max(NeuronValsRG(:,:,4),[],2))';
          
          %RespRG{i} = maxRG(logical(maxRG(:,:,2)),1,1)';

          RespRG{i} = (max(NeuronValsRG(respNeuronsMB<0.005,:,4),[],2))';
          pvalsRG{i} = respNeuronsRG(respNeuronsMB<0.005);

    else
        RespRG{i} = 0;
        pvalsRG{i} = 0;
       
    end


% 
%     p = NP.convertPhySorting2tIc(NP.recordingDir);
%     label = string(p.label');
%     goodU = p.ic(:,label == 'good');
%     goodUdepth = NP.chLayoutPositions(2,goodU(1,:));
%     uDepth = sin(deg2rad(data.Angle(ex))).*goodUdepth(goodNeurons);
%     verticalDepth = sin(deg2rad(data.Angle(ex))).*(data.Depth(ex)); %depth of unit along vertical axis
%     verticalDepthU= verticalDepth- uDepth;

%     OSIi = load(sprintf('Orientation-Tuning-Index-%s',NP.recordingName)).L';
%     DSIi = load(sprintf('Direction-Selectivity-Index-%s',NP.recordingName)).DSI';
% 
%     %theta{i} = load(sprintf('Angle-prefer-%s',NP.recordingName)).Theta;
%     OSI{i} = OSIi(goodNeurons);
%     DSI{i} = DSIi(goodNeurons);


    animalID{i} = data.Animal_ID(ex);

    tablePos{i} = ex;
    %DepthUnit{i} = Coor;
    i=i+1;

end
%%
%% Prepare data for summary plots
MB = cell2mat(RespMB);
RG = cell2mat(RespRG);

values = ([MB RG MB-RG])*(1000/300);

[uniqueStrings,~,indexAnimal] = unique([animalID{:}]);
color1 =[];
color2 =[];

j=1;

indexAnimal = sort(indexAnimal);
tablePosN = [];

for i = 1:length(animalID) %each color is an animal.
    
    color1 = [color1 zeros(1,length(RespMB{i}))+indexAnimal(i)];
    color2 = [color2 zeros(1,length(RespRG{i}))+indexAnimal(i)];
    tablePosN = [tablePosN repmat(tablePos{i},1,length(RespMB{i}))];
   
    j = j+1;


end

color = [color1 color2 color1];
tablePosN = [tablePosN tablePosN tablePosN];
NeuronID = [cell2mat(NeuronN) cell2mat(NeuronN) cell2mat(NeuronN)];


cats = categorical([ones(1,length(cell2mat(RespMB))) ones(1,length(cell2mat(RespRG)))+1 ones(1,length(cell2mat(RespMB)))+2]);



%%randomly select 10 neurons
MBh = length(MB)+length(RG)+find(MB-RG>0);

RGh = length(MB)+length(RG)+find(MB-RG<0);

randN = 5;
rng(42);
rMB = 1 + round((length(MBh) - 1) * rand(1, randN));
rng(42);
rRG = 1 + round((length(RGh) - 1) * rand(1, randN));

%selected neurons:

selecN = {[tablePosN((MBh(rMB)-length(MB)*2));NeuronID((MBh(rMB)-length(MB)*2))];...
    [tablePosN(RGh(rRG)-length(MB));NeuronID(RGh(rRG)-length(MB))]};

selecN{1}(2,1) = 218;
%% Plot responses

%randomNeur = 1 + round((length(cats) - 1) * rand(1, 10));

T = table(cats',values',color','VariableNames',{'cats','values','color'});

figure;

s = swarmchart(T.cats, T.values, 10, T.color,'filled','MarkerFaceAlpha',0.5);

hold on
swarmchart(T.cats([MBh(rMB(2))-length(MB)*2 MBh(rMB(2))-length(MB) MBh(rMB(2))]), ...
    T.values([MBh(rMB(2))-length(MB)*2 MBh(rMB(2))-length(MB) MBh(rMB(2))]), 50,...
    T.color([MBh(rMB(2))-length(MB)*2 MBh(rMB(2))-length(MB) MBh(rMB(2))]),...
    'filled','MarkerEdgeColor','r','LineWidth',2);

hold on
swarmchart(T.cats([RGh(rRG(5))-length(MB)*2 RGh(rRG(5))-length(MB) RGh(rRG(5))]),...
    T.values([RGh(rRG(5))-length(MB)*2 RGh(rRG(5))-length(MB) RGh(rRG(5))]), 50, ...
    T.color([RGh(rRG(5))-length(MB)*2 RGh(rRG(5))-length(MB) RGh(rRG(5))]),...
    'filled','MarkerEdgeColor','g','LineWidth',2);


myColormap = [
    1, 0, 0;    % Red
    0, 1, 0;    % Green
    0, 0, 1;    % Blue
    1, 1, 0;    % Yellow
    1, 0, 1;    % Magenta
    0, 1, 1;    % Cyan
    0, 0, 0;    % Black
]*0.7;
colormap(gca,myColormap)


% hold on;
% uniqueGroups = unique(T.color);  % Get unique group values
%  % Get the same colormap used in the plot
% 
% for i = [1 3]
%     % Plot invisible points for legend creation
%     scatter(nan, nan, 100, colors(i,:), 'filled');
%     
% end
% 
% 
% 
hold on
for i =1:3

    groupMean = mean(T.values(T.cats==categorical(i))); 

     plot([i-0.5, i+0.5], [groupMean, groupMean], 'r-', 'LineWidth', 2);  % Short horizontal line for the mean

end
% 
% legend({'PV35', 'PV103', 'SA5'}, 'Location', 'best');
grid on
set(gcf,'Color','w');%
xticklabels({'Moving','Static','Mov. - Stat.'})
ylabel('Response - Baseline (spikes/sec)')
ylim([-0.025 0.1])
cd('\\sil3\data\Large_scale_mapping_NP\lizards\SummaryFigs')
print(gcf,'ResponseMB-RG(red-High-MB)-(green-High-RG)-Color-Animal_1N)','-dpng')


%% Plot both signifficant responsive. 

MB = cell2mat(pvalsMB);
RG = cell2mat(pvalsRG);
values = [1-MB 1-RG];
color = [color1 color2];

cats = categorical([ones(1,length(cell2mat(RespMB))) ones(1,length(cell2mat(RespRG)))+1]);


T = table(cats',values',color','VariableNames',{'cats','values','color'});

figure;

s = swarmchart(T.cats, T.values, 10, T.color,'filled','MarkerFaceAlpha',0.5);

hold on
swarmchart(...
    T.cats([MBh(rMB(2))-length(MB)*2 MBh(rMB(2))-length(MB)]),...
    T.values([MBh(rMB(2))-length(MB)*2 MBh(rMB(2))-length(MB)]), 50,...
    T.color([MBh(rMB(2))-length(MB)*2 MBh(rMB(2))-length(MB)]),'filled','MarkerEdgeColor','r','LineWidth',2);

hold on
swarmchart(...
    T.cats([RGh(rRG(5))-length(MB)*2 RGh(rRG(5))-length(MB)]),...
    T.values([RGh(rRG(5))-length(MB)*2 RGh(rRG(5))-length(MB)]), 50, ...
    T.color([RGh(rRG(5))-length(MB)*2 RGh(rRG(5))-length(MB)]),...
    'filled','MarkerEdgeColor','g','LineWidth',2);


myColormap = [
    1, 0, 0;    % Red
    0, 1, 0;    % Green
    0, 0, 1;    % Blue
    1, 1, 0;    % Yellow
    1, 0, 1;    % Magenta
    0, 1, 1;    % Cyan
    0, 0, 0;    % Black
]*0.7;
colormap(gca,myColormap)


% hold on;
% uniqueGroups = unique(T.color);  % Get unique group values
%  % Get the same colormap used in the plot
% 
% for i = [1 3]
%     % Plot invisible points for legend creation
%     scatter(nan, nan, 100, colors(i,:), 'filled');
%     
% end
% 
% 
% 
hold on
for i =1:3

    groupMean = mean(T.values(T.cats==categorical(i))); 

     plot([i-0.5, i+0.5], [groupMean, groupMean], 'r-', 'LineWidth', 2);  % Short horizontal line for the mean

end
% 
% legend({'PV35', 'PV103', 'SA5'}, 'Location', 'best');
grid on
set(gcf,'Color','w');%
xticklabels({'Moving','Static'})
ylabel('1-(p-value)')
ylim([0.8 1.2])
cd('\\sil3\data\Large_scale_mapping_NP\lizards\SummaryFigs')
print(gcf,'ResponseMB-RG(red-High-MB)-(green-High-RG)-Color-Animal)','-dpng')


