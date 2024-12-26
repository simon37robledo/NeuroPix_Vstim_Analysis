%% Do the sumary graphs. Boxplot for moving ball. 

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

GoodRecordings =[1:18,40:43];%[1:20,28:32,40:48]; Perform bootstrap with SA5

%%
%GoodRecordings = 44;
animalN = 4;
i=1;
animalID = cell(1,animalN); %1
DepthUnit = cell(1,length(GoodRecordings)); %2
RespMB = cell(1,length(GoodRecordings));%3
RespRG = cell(1,length(GoodRecordings));%4
NeuronNMB = cell(1,length(GoodRecordings));%6
NeuronNRG =  cell(1,length(GoodRecordings));%6
tablePos = cell(1,length(GoodRecordings));%7
pvalsMB = cell(1,length(GoodRecordings));%8
pvalsRG = cell(1,length(GoodRecordings));%9
SpatialMB = cell(1,length(GoodRecordings));%10
SpatialRG = cell(1,length(GoodRecordings));%11
OSI =cell(1,length(GoodRecordings));%12
DSI = cell(1,length(GoodRecordings));%13
PreferAngle = cell(1,length(GoodRecordings));%14
EntropMB = cell(1,length(GoodRecordings));
EntropRG = cell(1,length(GoodRecordings));

zscoreRG= cell(1,length(GoodRecordings));%14
zscoreMB= cell(1,length(GoodRecordings));%14

MBsig =1; 
RGsig =0;
bothSig = 0;
allN = 0;
noEyeMoves = 0;
AllParams = cell(3,14);

close all

sign = 0.05;

N_bootstrap = 1000;

j = 1;

for ex = 40:43%GoodRecordings

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
%         MBsig =1; 
%         RGsig = 0;
        NeuronValsMB = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;
        respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        respNeuronsRG = load(sprintf('RectGrid-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        spaTuningMB = load(sprintf('Spatial-Tuning_index-MB-%s',NP.recordingName)).SpaTuning;
        tuningCurve = load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve;
        ZscoreNeuronsMB = load(sprintf('MovBall-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreU;
        ZscoreNeuronsRG = load(sprintf('RectGrid-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreU;
        

        if MBsig


            %sign = '0.005';
            zscoreMB{i} = ZscoreNeuronsMB(respNeuronsMB<sign);
            uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
            [preferDir dirInd] = max(tuningCurve,[],2);
            PreferAngle{i} = preferDir(respNeuronsMB<sign);
            RespMB{i} = max(NeuronValsMB(respNeuronsMB<sign,:,4),[],2)';
            NeuronNMB{i} = find(respNeuronsMB<sign);
            %RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
            pvalsMB{i} = respNeuronsMB(respNeuronsMB<sign);
            SpatialMB{i} = spaTuningMB(respNeuronsMB<sign);
            
            if noEyeMoves
                try
                    EntropMB{i} = load(sprintf('NEM-Entropies-MB-RF-respU-%s-%s.mat',sign,NP.recordingName)).entropies;
                catch
                    EntropMB{i} = [];
                end
            else
                try
                    EntropMB{i} = load(sprintf('Entropies-MB-RF-respU-%s-%s.mat',sign,NP.recordingName)).entropies;
                catch
                    EntropMB{i} = [];
                end

            end

            %DepthUnit{i} =
            colors = repmat([0 0 0],size(goodU,2),1);
            if sum(respNeuronsMB<0.005) >0
                if j ==1
                    [Coor, Fig] = neuronLocation(NP,data(ex,:),goodU(:,respNeuronsMB<sign),1,0,colors, [], 1);

                else

                    [Coor,~] = neuronLocation(NP,data(ex,:),goodU(:,respNeuronsMB<sign),1,0,colors, Fig, 1);

                end
                j =j+1;

            end



      

        end

        if RGsig

            zscoreMB{i} = ZscoreNeuronsMB(respNeuronsRG<sign);
            uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
            [preferDir dirInd] = max(tuningCurve,[],2);
            PreferAngle{i} = preferDir;
            RespMB{i} = max(NeuronValsMB(respNeuronsRG<sign,:,4),[],2)';
            NeuronNMB{i} = find(respNeuronsRG<sign);
            %RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
            pvalsMB{i} = respNeuronsRG(respNeuronsRG<sign);
            SpatialMB{i} = spaTuningMB(respNeuronsRG<sign);

        end
        if allN
            
            uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
            [preferDir dirInd] = max(tuningCurve,[],2);
            PreferAngle{i} = preferDir;
            RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
            NeuronNMB{i} = 1:length(respNeuronsRG);
            %RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
            pvalsMB{i} = respNeuronsRG;
            SpatialMB{i} = spaTuningMB;
            zscoreMB{i} =0;

        end

        if bothSig

            zscoreMB{i} = ZscoreNeuronsMB(respNeuronsRG<sign & respNeuronsMB<sign);
            uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
            [preferDir dirInd] = max(tuningCurve,[],2);
            PreferAngle{i} = preferDir;
            RespMB{i} = max(NeuronValsMB(respNeuronsRG<sign & respNeuronsMB<sign,:,4),[],2)';
            NeuronNMB{i} = find(respNeuronsRG<0.005 & respNeuronsMB<0.005);
            %RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
            pvalsMB{i} = respNeuronsRG(respNeuronsRG<sign & respNeuronsMB<sign);
            SpatialMB{i} = spaTuningMB(respNeuronsRG<sign & respNeuronsMB<sign);


        end

    else
        goodNeuronsMB{i} =0;
        PreferAngle{i} = 0;
        RespMB{i} = 0;
        NeuronNMB{i} =0;
        pvalsMB{i} =0;
        SpatialMB{i} =0;
        EntropMB{i} = 0;


    end


    if contains(data.VS_ordered(ex),"RG") &&  contains(data.VS_ordered(ex),"MB")
%         RGsig =1;
%         MBsig = 0;
        NeuronValsRG= load(sprintf('RectGrid-NeuronRespCat-%s',NP.recordingName)).NeuronVals;
        spaTuningRG = load(sprintf('Spatial-Tuning_index-RG-%s',NP.recordingName)).SpaTuning;
        respNeuronsRG = load(sprintf('RectGrid-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        ZscoreNeuronsRG = load(sprintf('RectGrid-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreU;
        
        if MBsig
            
            RespRG{i} = (max(NeuronValsRG(respNeuronsMB<sign,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG(respNeuronsMB<sign);
            SpatialRG{i} = spaTuningRG(respNeuronsMB<sign);
            NeuronNRG{i} = find(respNeuronsMB<sign);
            zscoreRG{i} = ZscoreNeuronsRG(respNeuronsMB<sign);

            if noEyeMoves
                try
                    EntropRG{i} = load(sprintf('NEM-Entropies-RG-RF-respU-%s.mat',NP.recordingName)).entropies;
                catch
                    EntropRG{i} = [];
                end
            else
                try
                    EntropRG{i} = load(sprintf('Entropies-RG-RF-respU-%s-%s.mat',NP.recordingName)).entropies;
                catch
                    EntropRG{i} = [];
                end
                
            end

            EntropRG{i} = entropiesRG(respNeuronsMB<sign);
        end


        if RGsig
            RespRG{i} = (max(NeuronValsRG(respNeuronsRG<sign,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG(respNeuronsRG<sign);
            SpatialRG{i} = spaTuningRG(respNeuronsRG<sign);
            NeuronNRG{i} = find(respNeuronsRG<sign);
            zscoreRG{i} = ZscoreNeuronsRG(respNeuronsRG<sign);
        end


        if allN 
            RespRG{i} = (max(NeuronValsRG(:,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG;
            SpatialRG{i} = spaTuningRG;
            NeuronNRG{i} = 1:length(respNeuronsRG);
            zscoreRG{i} =0;
        end

        if bothSig
            RespRG{i} = (max(NeuronValsRG(respNeuronsRG<sign & respNeuronsMB<sign,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG(respNeuronsRG <sign & respNeuronsMB<sign);
            SpatialRG{i} = spaTuningRG(respNeuronsRG<sign & respNeuronsMB<sign);
            NeuronNRG{i} = find(respNeuronsRG<sign & respNeuronsMB<sign);
            zscoreRG{i} = ZscoreNeuronsRG(respNeuronsRG<sign & respNeuronsMB<sign);
             EntropRG{i} = entropiesRG(respNeuronsRG<sign & respNeuronsMB<sign);
        end

    else
        RespRG{i} = 0;
        pvalsRG{i} = 0;
        SpatialRG{i} = 0;
        NeuronNRG{i} = 0;
        EntropRG{i} = 0;

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
%
%% Prepare data for summary plots
MB = cell2mat(RespMB);
RG = cell2mat(RespRG);

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
NeuronID = [cell2mat(NeuronNMB) cell2mat(NeuronNRG) cell2mat(NeuronNMB)];


cats = categorical([ones(1,length(cell2mat(RespMB))) ones(1,length(cell2mat(RespRG)))+1 ones(1,length(cell2mat(RespMB)))+2]);



%% randomly select 10 neurons
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

%selecN{1}(2,1) = 218;
%% Plot responses
MB = cell2mat(RespMB);
RG = cell2mat(RespRG);

values = ([MB RG MB-RG])*(1000/300);
%randomNeur = 1 + round((length(cats) - 1) * rand(1, 10));

T = table(cats',values',color','VariableNames',{'cats','values','color'});

figure;

s = swarmchart(T.cats, T.values, 10, T.color,'filled','MarkerFaceAlpha',0.5);

hold on %%Chose before, 2 and 5;
chooseNb = 1:2;
swarmchart(T.cats([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]), ...
    T.values([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]), 50,...
    T.color([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]),...
    'filled','MarkerEdgeColor','r','LineWidth',2);

hold on
chooseNr = 1:2;

swarmchart(T.cats([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]),...
    T.values([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]), 50, ...
    T.color([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]),...
    'filled','MarkerEdgeColor','g','LineWidth',2);
% 
% 
% myColormap = [
%     1, 0, 0;    % Red
%     0, 1, 0;    % Green
%     0, 0, 1;    % Blue
%     1, 1, 0;    % Yellow
%     1, 0, 1;    % Magenta
%     0, 1, 1;    % Cyan
%     0, 0, 0;    % Black
% ]*0.7;
% colormap(gca,myColormap)


myColormap = [
    1, 0, 0;    % Red
    0, 1, 0;    % Green
    0, 0, 1;    % Blue
    1, 0, 1;    % Magenta
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
%ylim([-0.025 0.1])
cd('\\sil3\data\Large_scale_mapping_NP\lizards\SummaryFigs')
print(gcf,'SIg(MB)-ResponseMB-RG(red-High-MB)-(green-High-RG)-Color-Animal_1N)','-dpng')


%% Plot Z-scores
MB = cell2mat(zscoreMB);
RG = cell2mat(zscoreRG);

eMB = cell2mat(EntropMB);
eRG = cell2mat(EntropRG);

values = ([MB RG MB-RG]);

values = ([eMB eRG eMB-eRG]);

%%randomly select 10 neurons but from animals that have full receptive
%%field tuning == PV35 & PV139 19:20 40:43
%goodRFr = cell2mat(zscoreMB(1))-cell2mat(zscoreRG(19:24));
goodRFr = cell2mat(EntropMB(1))-cell2mat(EntropRG(1));
%MBh = length(MB)+length(RG)+length(cell2mat(zscoreMB(1)))+find(goodRFr<0);
MBh = length(MB)+length(RG)+find(goodRFr<0);
RGh = length(MB)+length(RG)+find(goodRFr>0);

seed1 = 42;

seed2 = 43;
randN = 5;
rng(42);
rMB = 1 + round((length(MBh) - 1) * rand(1, randN));
rng(42);
rRG = 1 + round((length(RGh) - 1) * rand(1, randN));

%selected neurons:

%selecN = {[tablePosN((MBh(rMB)-length(MB)*2));NeuronID((MBh(rMB)-length(MB)*2))];...
  %  [tablePosN(RGh(rRG)-length(MB));NeuronID(RGh(rRG)-length(MB))]};


%randomNeur = 1 + round((length(cats) - 1) * rand(1, 10));

T = table(cats',values',color','VariableNames',{'cats','values','color'});

figure;

s = swarmchart(T.cats, T.values, 10, T.color,'filled','MarkerFaceAlpha',0.5);

hold on %%Chose before, 2 and 5;
chooseNb = 1;
swarmchart(T.cats([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]), ...
    T.values([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]), 70,...
    T.color([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]),...
    'filled','MarkerEdgeColor','r','LineWidth',2);

% hold on
% chooseNr = 1;
% swarmchart(T.cats([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]),...
%     T.values([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]), 50, ...
%     T.color([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]),...
%     'filled','MarkerEdgeColor','g','LineWidth',2);
% 
% 
% myColormap = [
%     1, 0, 0;    % Red
%     0, 1, 0;    % Green
%     0, 0, 1;    % Blue
%     1, 1, 0;    % Yellow
%     1, 0, 1;    % Magenta
%     0, 1, 1;    % Cyan
%     0, 0, 0;    % Black
% ]*0.7;
% colormap(gca,myColormap)


% myColormap = [
%     1, 0, 0;    % Red
%     0, 1, 0;    % Green
%     0, 0, 1;    % Blue
%     1, 0, 1;    % Magenta
%     0, 0, 0;    % Black
% ]*0.7;
% colormap(gca,myColormap)

% myColormap = [
%     0, 0, 1;    % Blue
%     1, 0, 1;    % Magenta
%     0, 0, 0;    % Black
% ]*0.7;
% colormap(gca,myColormap)



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
ylabel('Entropy')
yL = ylim;
%ylim([yL(1) 60])
cd('\\sil3\data\Large_scale_mapping_NP\lizards\SummaryFigs')
%c = colorbar;
% title(c, 'Animals (N=5)');%caxis([1 5])
%c.Ticks = [linspace(1.5,4.5,5)];
% c.Ticks = [linspace(1.5,2.5,animalN)];
% c.TickLabels = {'PV67','PV139',  'PV35'};

print(gcf,'NEM-AWAKESIg(MB)-entropy-MB-RG(red-High-MB)-(green-High-RG)-Color-Animal_1N)','-dpng')
% 
 %% Plot both signifficant responsive. 

MB = cell2mat();
RG = cell2mat(pvalsRG);
values = [1-MB 1-RG];

cats = categorical([ones(1,length(cell2mat(RespMB))) ones(1,length(cell2mat(RespRG)))+1]);

color1 =[];
color2 =[];
for i = 1:length(animalID) %each color is an animal.
    
    color1 = [color1 zeros(1,length(RespMB{i}))+indexAnimal(i)];
    color2 = [color2 zeros(1,length(RespRG{i}))+indexAnimal(i)];
    tablePosN = [tablePosN repmat(tablePos{i},1,length(RespMB{i}))];
   
    j = j+1;


end

color = [color1 color2 color1];


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

%% Plot spatial tuning
MB = cell2mat(SpatialMB);
RG = cell2mat(SpatialRG);
%RG = RG(cell2mat(pvalsRG)<0.005);

values = ([MB RG]);
cats = categorical([ones(1,length(MB)) ones(1,length(RG))+1]);
%randomNeur = 1 + round((length(cats) - 1) * rand(1, 10));

color1 = [];
color2 = [];
tablePosN1 = [];
tablePosN2 = [];
for i = 1:length(animalID) %each color is an animal.
    
    color1 = [color1 zeros(1,length(RespMB{i}))+indexAnimal(i)];
    color2 = [color2 zeros(1,length(RespRG{i}))+indexAnimal(i)];
    tablePosN1 = [tablePosN1 repmat(tablePos{i},1,length(RespMB{i}))];
    tablePosN2 = [tablePosN2 repmat(tablePos{i},1,length(RespRG{i}))];
    j = j+1;


end

color = [color1 color2];

tablePosN = [tablePosN1 tablePosN2];

NeuronID = {cell2mat(NeuronNMB) cell2mat(NeuronNRG)};



T = table(cats',values',color','VariableNames',{'cats','values','color'});

figure;

s = swarmchart(T.cats, T.values, 10, T.color,'filled','MarkerFaceAlpha',0.5);



%%randomly select 10 neurons
SPIhMB = find(MB>0.8);
SPIhRG = find(RG>0.8);
SPIlMB = find(MB<0.5);
SPIlRG = find(RG<0.5);
randN = 1;
rng(42);
rHMB = 1 + round((length(SPIhMB) - 1) * rand(1, randN));
rHRG = 1 + round((length(SPIhRG) - 1) * rand(1, randN));
rng(42);
rLMB = 1 + round((length(SPIlMB) - 1) * rand(1, randN));
rLRG = 1 + round((length(SPIlRG) - 1) * rand(1, randN));
%selected neurons:


%%higher SPT = first cell (row 1 MB, row 2 RG) / Lower SPT second cell (row 1 MB, row 2 RG)
% selecN = {[tablePosN1(SPIhMB(rHMB)) NeuronID{1}(SPIhMB(rHMB));tablePosN2(SPIhRG(rHRG)) NeuronID{2}(SPIhRG(rHRG))];...
%     [tablePosN1(SPIlMB(rLMB)) NeuronID{1}(SPIlMB(rLMB));tablePosN2(SPIlRG(rLRG)) NeuronID{2}(SPIlRG(rLRG))]};

selecN = {[tablePosN1(SPIhMB(rHMB)) NeuronID{1}(SPIhMB(rHMB))];...
    [tablePosN1(SPIlMB(rLMB)) NeuronID{1}(SPIlMB(rLMB))]};

%%Plot selected neurons

hold on %%Chose before, 2 and 5;
chooseNb = 1;

%Table representation of selected neurons for the two cathegories
%TRH = [SPIhMB(rHMB(chooseNb)) SPIhRG(rHRG(chooseNb))+length(MB)];
TRH = [SPIhMB(rHMB(chooseNb)) SPIhMB(rHMB(chooseNb))+length(MB)];
swarmchart(T.cats(TRH), ...
    T.values(TRH), 50,...
    T.color(TRH),...
    'filled','MarkerEdgeColor','r','LineWidth',2);

hold on
chooseNr = 1;
%TRL = [SPIlMB(rLMB((chooseNr))) SPIlRG(rLRG(chooseNr))+length(MB)];
TRL = [SPIlMB(rLMB((chooseNr))) SPIlMB(rLMB((chooseNr)))+length(MB)];
swarmchart(T.cats(TRL),...
    T.values(TRL), 50, ...
    T.color(TRL),...
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
for i =1:length(unique(T.cats))

    groupMean = mean(T.values(T.cats==categorical(i))); 

     plot([i-0.5, i+0.5], [groupMean, groupMean], 'r-', 'LineWidth', 2);  % Short horizontal line for the mean

end
% 
% legend({'PV35', 'PV103', 'SA5'}, 'Location', 'best');
grid on
set(gcf,'Color','w');%
xticklabels({'Moving','Static','Mov. - Stat.'})
ylabel('Spatial Tuning Index')
%ylim([-0.025 0.1])
cd('\\sil3\data\Large_scale_mapping_NP\lizards\SummaryFigs')
print(gcf,'SIg(MB)-SpatialTuning-Color-Animal_1N)','-dpng')

