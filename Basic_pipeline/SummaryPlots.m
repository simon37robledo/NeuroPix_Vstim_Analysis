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

N_bootstrap = 1000;

for ex = GoodRecordings

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


    if contains(data.VS_ordered(ex),"MB")
        NeuronValsMB = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;

        respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        goodNeuronsMB{i} = find(respNeuronsMB<0.005);
        tuningCurve = load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve;
        uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
        [preferDir dirInd] = max(tuningCurve,[],2);
        PreferAngle{i} = preferDir(respNeuronsMB<0.005);
        RespMB{i} = max(NeuronValsMB(respNeuronsMB<0.005,:,4),[],2)';

    else
        goodNeuronsMB{i} =0;
        PreferAngle{i} = 0;
        RespMB{i} = 0;

    end


    if contains(data.VS_ordered(ex),"RG") &&  contains(data.VS_ordered(ex),"MB")

          NeuronValsRG= load(sprintf('RectGrid-NeuronRespCat-%s',NP.recordingName)).NeuronVals;
          %EnoughTrials = NeuronValsRG(respNeuronsMB<0.005,:,[4 7]);

          maxRG = (max(NeuronValsRG(respNeuronsMB<0.005,:,[4 7]),[],2));
          
          RespRG{i} = maxRG(logical(maxRG(:,:,2)),1,1)';
    else
        RespRG{i} = 0;
       
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
    %DepthUnit{i} = Coor;
    i=i+1;

end
%%
%% Swarm plots 

values = ([RespMB{:} RespRG{:}])*(1000/duration);

[uniqueStrings,~,indexAnimal] = unique([animalID{:}]);
color1 =[];
color2 =[];

j=1;

indexAnimal = sort(indexAnimal);

for i = 1:length(animalID) %each color is an animal.
    
    color1 = [color1 zeros(1,length(RespMB{i}))+indexAnimal(i)];
    color2 = [color2 zeros(1,length(RespRG{i}))+indexAnimal(i)];
    j = j+1;

end

color = [color1 color2];

cats = categorical([ones(1,length(cell2mat(RespMB))) ones(1,length(cell2mat(RespRG)))+1]);

xticklabels({'MB','RG'})


T = table(cats',values',color','VariableNames',{'cats','values','color'});

figure;

swarmchart(T.cats, T.values, 10, T.color,'filled');

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
for i =1:2

    groupMean = mean(T.values(T.cats==categorical(i))); 

     plot([i-0.5, i+0.5], [groupMean, groupMean], 'r-', 'LineWidth', 2);  % Short horizontal line for the mean

end
% 
% legend({'PV35', 'PV103', 'SA5'}, 'Location', 'best');
grid on
set(gcf,'Color','w');%
xticklabels({'Moving Ball','Rectangle Grid'})
ylabel('Spikes/sec - baseline spkRate')
ylim([-0.01 0.1])