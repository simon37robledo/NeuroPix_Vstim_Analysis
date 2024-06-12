%%% Rect Grid Analysis

path = '\\sil3\data\\Large_scale_mapping_NP\lizards\PV139\PV139_Experiment_6_2_24\Insertion1\catgt_PV139_Experiment_6_2_24_1_g0';


basic_pathPV102 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV102';
expPV102 = 'PV102_experiment_18_7_23';

basic_pathPV103 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV103';
expPV103 = 'PV103_Experiment_12_6_23';

basic_pathPV67 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV67';
expPV67= 'PV67_experiment_5_7_23';

basic_pathPV27 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV27';
expPVPV27 = 'PV27_Experiment_25_6_23';

basic_pathPV139 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV139';
expPV139 = 'PV139_Experiment_6_2_24';

basic_pathPV59 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV59';
expPV59 = 'PV59_Experiment_20_2_24';

basic_pathPV32 = '\\sil3\data\Large_scale_mapping_NP\lizards\PV32';
expPV32 = 'PV32_Experiment_18_3_24';

basicPathA = {basic_pathPV67,basic_pathPV103,basic_pathPV27,basic_pathPV139,basic_pathPV59,basic_pathPV32};

expA = {expPV67,expPV103,expPVPV27,expPV139,expPV59,expPV32};

RGttl_Index = {3,1,1,[2,1],2,2};

a =1;
in ="4";

path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
                in+string(filesep)+"catgt_"+string(expA{a})+"_"+in+"_g0");

NP = NPAPRecording(path); 
cd(path)

newExp =false;

%%


%%%Moving ball inputs, NP, stimOn, stimOff, plots, neuron, ttlIndex


patternIndex = strfind(string(NP.recordingDir), "\catgt");

endIndex = patternIndex(1)-1;
stimDir = string(NP.recordingDir);
stimDir = extractBetween(stimDir,1,endIndex);

file = dir (stimDir);
filenames = {file.name};

file = dir (stimDir);
filenames = {file.name};
rectFiles = filenames(contains(filenames,"rectGrid"));
positionsMatrix = [];
offsets = [];
sizes = [];
seqMatrix = [];

j =1;
if size(rectFiles) ~= [0 0]

    for i = rectFiles
        rect= load(stimDir+"\"+string(i));

        if newExp 
            %%New exp
            seqMatrix = [seqMatrix cell2mat(rect.VSMetaData.allPropVal(18))];
            sizes = [sizes cell2mat(rect.VSMetaData.allPropVal(23))];
            stimDurStats = cell2mat(rect.VSMetaData.allPropVal(42))*1000;
            interStimStats = cell2mat(rect.VSMetaData.allPropVal(32))*1000;
        else
            %Old exp
            seqMatrix = [seqMatrix cell2mat(rect.VSMetaData.allPropVal(16))];
           
            stimDurStats = cell2mat(rect.VSMetaData.allPropVal(38))*1000;
            interStimStats = cell2mat(rect.VSMetaData.allPropVal(28))*1000;
        end

        j = j+1;
    end
    disp('Visual stats extracted!')
else
    disp('Directory does not exist!');
end

if newExp
    positionsMatrix = [cell2mat(rect.VSMetaData.allPropVal(16)),cell2mat(rect.VSMetaData.allPropVal(17))];%NewExp
    rectSizes = cell2mat(rect.VSMetaData.allPropVal(22));
else
    positionsMatrix = [cell2mat(rect.VSMetaData.allPropVal(14)),cell2mat(rect.VSMetaData.allPropVal(15))]; %OldExp
    rectSizes = cell2mat(rect.VSMetaData.allPropVal(22));
    sizes = repmat(cell2mat(rect.VSMetaData.allPropVal(5)),1,length(seqMatrix));
end

%% Triggers:

ttlInd = RGttl_Index{a};
[stimOn stimOff] = NPdiodeExtract(NP,0,"rectGrid",ttlInd,16,0); 

[stimOn stimOff] = NPdiodeExtract(NP,0,"rectGrid",ttlInd,16,0); 
missedT = find((stimOff-stimOn)<500); %Missed tri als

stimDur = mean(stimOff'-stimOn');




A = [stimOn seqMatrix' sizes'];

[C indRG]= sortrows(A,[2 3]);

%Sort directions:

directimesSorted = C(:,1)';


%%

cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

%Good units

GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

p = NP.convertPhySorting2tIc(NP.recordingDir);

%Select good units
label = string(p.label');

goodU = p.ic(:,label == 'good');


%%
% Build raster plots per unit

cd(NP.recordingDir)
if ~exist(path+"\Figs",'dir') 
    mkdir Figs
end
cd(NP.recordingDir + "\Figs")

bin=20;
win=stimDur+stimInter;
preBase = round(stimInter/20)*10;

[Mr]=BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(win/bin));

[nT,nN,nB] = size(Mr);
%indRG --> sorted infexes 

trialsPerCath = length(seqMatrix)/(length(unique(seqMatrix)));

PositionsTotal = positionsMatrix(seqMatrix(indRG),:);


[posS,indexX] = sortrows(PositionsTotal,1); %Sort first dimension because tile layout moves through columns

MrS = Mr(indexX,:,:);

nSize = length(unique(sizes));



for u =1:length(goodU)
    t = tiledlayout(sqrt(max(seqMatrix)), sqrt(max(seqMatrix)),'TileSpacing','tight');

    j=1;

    if trialsPerCath>20
        mergeTrials = 5;

        Mr2 = zeros(nT/mergeTrials,nB);

        for i = 1:mergeTrials:nT

            meanb = mean(squeeze(MrS(i:min(i+mergeTrials-1, end),u,:)),1);

            Mr2(j,:) = meanb;

            j = j+1;

        end
    else
        Mr2=MrS(:,u,:);
        mergeTrials =1;
    end


    [T,B] = size(Mr2);

    for i = 1:trialsPerCath/mergeTrials:T
        %Build raster
        M = Mr2(i:min(i+trialsPerCath/mergeTrials-1, end),:);
        [nTrials,nTimes]=size(M);
        nexttile
        imagesc((1:nTimes)*bin,1:nTrials,squeeze(M));colormap(flipud(gray(64)));
        xline(preBase, '--', LineWidth=2, Color="#77AC30");
        xline((stimDur+preBase), '--', LineWidth=2, Color="#0072BD");
        xticks([preBase round(round(stimDur/100))*100+preBase]);
        if nSize >1
            yline([trialsPerCath/mergeTrials/nSize:trialsPerCath/mergeTrials/nSize:trialsPerCath/mergeTrials-1]+0.5,LineWidth=1)
        end
        caxis([0 1]);
        set(gca,'YTickLabel',[]);

        if i < T - (trialsPerCath/mergeTrials)*max(positionsMatrix(:))-1
            set(gca,'XTickLabel',[]);

        end
    end
    fig = gcf;
    set(fig, 'Color', 'w');
    colorbar
    % Set the color of the figure and axes to black
    title(t,sprintf('Rect-GRid-raster-U%d',u))
    ylabel(t,sprintf('%d trials',nTrials*mergeTrials))
    fig.Position = [227         191        1413         781];
    %prettify_plot
    print(fig,sprintf('%s-rect-GRid-raster-U%d.png',NP.recordingName,u),'-dpng')
    close
end

%% Build heat map

cd(NP.recordingDir)
if ~exist(path+"\Figs",'dir')
    mkdir Figs
end
cd(NP.recordingDir + "\Figs")

trialDiv  = length(seqMatrix)/length(unique(seqMatrix))/nSize;

offsetR=400; %ms for offset response

bin =10;
[Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round(stimDur+offsetR)/bin);

[nT,nN,NB] = size(Mr);

MrC = zeros(round(nT/trialDiv),nN, floor(nB*2.5/4)-floor(nB/10)+floor(offsetR/bin));


%%Create summary of identical trials

for u = 1:length(goodU)
    j=1;

    for i = 1:trialDiv:nT

        meanRon = mean(squeeze(Mr(i:i+trialDiv-1,u,floor(nB/10)+1:floor(nB*2.5/4))));

        meanRoff =  mean(squeeze(Mr(i:i+trialDiv-1,u,floor(stimDur/bin)+1:end)));

        MrC(j,u,:) = [meanRon meanRoff]; %Combine on and off response

        j = j+1;

    end
end

MrMean = mean(MrC,3);


[Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin));

Nb2=  mean(Mb,3);

Nbase = mean(Mb,[1 3]);


%figure;imagesc(squeeze(MrC(:,18,:)));

points = {'.','+','*','o'};
%%
cd(NP.recordingDir + "\Figs")
respU = [138,134,127,124,123,121,118,112]; %PV67-1 %check Offset response also. 
respU = [51,28];
% Create a meshgrid of coordinates

screenSide = rect.VSMetaData.allPropVal{64}; %Same as moving ball
[x, y] = meshgrid(1:screenSide, 1:screenSide);
%screenSide = (max(rect.VSMetaData.allPropVal{21,1}.Y4{1,4},[],'all'))

%exU =[1];

cd

RFuT = zeros(screenSide,screenSide,length(u));

for u = respU

    j=1;

   %u =12;

  %  RFu = zeros(screenSide,screenSide,length(C)/);

  pxyScreen = zeros(screenSide,screenSide);

    VideoScreen = zeros(screenSide,screenSide,length(C)/trialDiv);

    if newExp
        prop = 21;
    else
        prop =19;
    end
   
    for i = 1:trialDiv:length(C)

        xyScreen = zeros(screenSide,screenSide); %%Make calculations if sizes>1 and if experiment is new and the shape is a circle. 
       
        if nSize>1 && string(rect.VSMetaData.allPropVal{7}) == "circle"
            Xc = round((rect.VSMetaData.allPropVal{prop,1}.X2{1,C(i,3)}(C(i,2))-rect.VSMetaData.allPropVal{prop,1}.X1{1,C(i,3)}(C(i,2)))/2+rect.VSMetaData.allPropVal{prop,1}.X1{1,C(i,3)}(C(i,2)));%...
            %-min(rect.VSMetaData.allPropVal{21,1}.Y2{1,4},[],'all'))+conversion;
            Yc = round((rect.VSMetaData.allPropVal{prop,1}.Y4{1,C(i,3)}(C(i,2))-rect.VSMetaData.allPropVal{prop,1}.Y1{1,C(i,3)}(C(i,2)))/2+rect.VSMetaData.allPropVal{prop,1}.Y1{1,C(i,3)}(C(i,2)));%...
            %-min(rect.VSMetaData.allPropVal{21,1}.Y2{1,4},[],'all'))+conversion;

            r = round((rect.VSMetaData.allPropVal{prop,1}.X2{1,C(i,3)}(C(i,2))-rect.VSMetaData.allPropVal{prop,1}.X1{1,C(i,3)}(C(i,2)))/2);

            % maxX = rect.VSMetaData.allPropVal{21,1}.X2{1,1};


            % Calculate the distance of each point from the center
            distances = sqrt((x - Xc).^2 + (y - Yc).^2);

            % Set the values inside the circle to 1 (or any other value you prefer)
            xyScreen(distances <= r) = 1;

            %
            %                 hold on; plot(rect.VSMetaData.allPropVal{21,1}.X1{1,C(i,3)}(C(i,2)),rect.VSMetaData.allPropVal{21,1}.Y1{1,C(i,3)}(C(i,2)),points{C(i,3)},MarkerSize=10)%,...
            %                 hold on; plot(rect.VSMetaData.allPropVal{21,1}.X2{1,C(i,3)}(C(i,2)),rect.VSMetaData.allPropVal{21,1}.Y2{1,C(i,3)}(C(i,2)),points{C(i,3)},MarkerSize=10)
            %                 hold on; plot(rect.VSMetaData.allPropVal{21,1}.X3{1,C(i,3)}(C(i,2)),rect.VSMetaData.allPropVal{21,1}.Y3{1,C(i,3)}(C(i,2)),points{C(i,3)},MarkerSize=10)%,...
            %                 hold on; plot(rect.VSMetaData.allPropVal{21,1}.X4{1,C(i,3)}(C(i,2)),rect.VSMetaData.allPropVal{21,1}.Y4{1,C(i,3)}(C(i,2)),points{C(i,3)},MarkerSize=10)%,...
            %                 hold on; plot(Xc,Yc,points{C(i,3)},MarkerSize=10);

            VideoScreen(:,:,j) = xyScreen;

        elseif nSize ==1 %%Add the possibility of several sizes but square. 

            X2 = round(rect.VSMetaData.allPropVal{prop,1}.X2(C(i,2)));
            X1 = round(rect.VSMetaData.allPropVal{prop,1}.X1(C(i,2)));
            X3 = round(rect.VSMetaData.allPropVal{prop,1}.X3(C(i,2)));
            X4 = round(rect.VSMetaData.allPropVal{prop,1}.X4(C(i,2)));
            X = [X1,X2,X3,X4];

            Y4 = round(rect.VSMetaData.allPropVal{prop,1}.Y4(C(i,2)));
            Y1 = round(rect.VSMetaData.allPropVal{prop,1}.Y1(C(i,2)));
            Y3 = round(rect.VSMetaData.allPropVal{prop,1}.Y3(C(i,2)));
            Y2 = round(rect.VSMetaData.allPropVal{prop,1}.Y2(C(i,2)));
            Y = [Y1,Y2,Y3,Y4];

            mask = poly2mask(X, Y, screenSide, screenSide);
            xyScreen(mask) =1;
            %figure;imagesc(xyScreen)

            VideoScreen(:,:,j) = xyScreen;

        end
        % x = zeros(1800,1800,100*length(C)/trialDiv);

       pxyScreen = pxyScreen+xyScreen;

       %[Mc]=convn(squeeze(M)',fspecial('gaussian',[1 5],3),'same');

       %M(1,1,:) = Mc;

       %conV = convn(VideoScreen,M,'same');

       %[MaxIm ind] = max(conV(Xc,Yc,:)); %%Delay according to max

       %conVM = conV(:,:,ind);

       %RFu = RFu+conVM;

       j=j+1;
        
%        figure; plot(1:length(squeeze(M)),squeeze(M)')
% 
%        figure; plot(1:length(squeeze(M)),squeeze(Mc)')
    end

  %  M = MrMean(:,u)'./Nbase(u);

    epsilon = 0.01;

    denom = mad(Nb2(u),0)+epsilon; %mean(Mb,0)+epsilon; %

   % mSpk = mean(spkRateR);

    M = ( MrMean(:,u)' - (Nbase(u) + MrMean(:,u)')./2)./denom;


    RFu = mean(bsxfun(@times, VideoScreen, reshape(M, 1, 1, [])),3);
    %figure;imagesc(RFu)
    fig = figure;imagesc(RFu./pxyScreen)
    caxis([0 0.02]);
    colorbar; max(RFu./pxyScreen,[],'all')
    xlabel('X pixels')
    ylabel('Y pixels')
    title(sprintf('RFu-%d-Static',u))
    prettify_plot
    
    print(fig,sprintf('%s-RFu-%d-Static.png',NP.recordingName,u),'-dpng')
    close
  

end












