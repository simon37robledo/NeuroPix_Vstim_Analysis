%% Do the sumary graphs. Boxplot for moving ball. 

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

allGoodRec = [1:21, 28:36, 40:55];
GoodRecordings =[8:21,40:43, 49:54];%Anesthetized
FFFrecordingsA = [15:18,40:43]; %anesthetized
SDGrecordingsA = [8:14,40:43,49:54]; %anesthetized 
%%% MB analysis PV97 done
%%% SDG analysis PV97 done
%%% RG analysis PV97 done
%%% FFF analysis PV97 done

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
unitCoorX = cell(1,length(GoodRecordings)); %medial-lateral
unitCoorY = cell(1,length(GoodRecordings)); %frontal-caudal
unitCoorZ = cell(1,length(GoodRecordings)); %superior-inferior
animalName ={};

%Neuron counts
nTFFF = 0;
nTSDGs = 0;
nTSDGm = 0;
nTMB = 0;
nTRG = 0;
totalN =0;

MBtotal = 0;
RGtotal =0;
SDGmtotal =0;
SDGstotal = 0;
FFFtotal =0;

onlyMB = 0;
onlyRG = 0;
onlySDGm = 0;
onlySDGs = 0;
SDGmAndRG = 0;
MBandRG = 0;

SDGmAndSDGs = 0;
SDGmAndMB = 0;
MBandSDGs = 0;
MBandSDGmsndSDGs = 0;
MB = [];
RG =[];

BombcelledNeurons = 1;


%%%What neurons to select? Only one active at the same time
MBsig =1; %%Neurons responsive to MB
MBsigAnd =0; %%Neurons responsive to MB and each other stimuli
RGsig =0;
bothSig = 0;
allN = 0;
eachSign = 0;
noEyeMoves = 0;
AllParams = cell(3,14);

%%% Initialize z-scores accordying to venn diagram
ZSonlyMB = []; 
ZSonlySDGm = [];
ZSonlyRG = [];
ZSonlySDGs = [];

%and color vectors:
MBc =[];
SDGmC = [];
SDGsC = [];
RGc =[];

%%%What to plot?
BrainPlots=0;
DSIplot =0;
Zscoreplot = 0;
Venn = 0;

close all

sign = 0.005;

N_bootstrap = 1000;

%
j = 1;
animal =0;

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

%     if BombcelledNeurons
% 
%         somaticNeur 
% 
%     end


     %%%%%%%%%%%%%%%%%%% MB 

    if contains(data.VS_ordered(ex),"MB")
%         MBsig =1; 
%         RGsig = 0;
        NeuronValsMB = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;
        respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        totalN = totalN+length(respNeuronsMB);
        respNeuronsRG = load(sprintf('RectGrid-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        %spaTuningMB = load(sprintf('Spatial-Tuning_index-MB-%s',NP.recordingName)).SpaTuning;
        tuningCurve = load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve;
        ZscoreNeuronsMB = load(sprintf('MovBall-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreU;
        ZscoreNeuronsRG = load(sprintf('RectGrid-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreU;
        

        if MBsig || MBsigAnd


            %sign = '0.005';
            zscoreMB{i} = ZscoreNeuronsMB(respNeuronsMB<=sign);
            nTMB = nTMB+length(ZscoreNeuronsMB(respNeuronsMB<=sign));
            uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
            [preferDir dirInd] = max(tuningCurve,[],2);
            PreferAngle{i} = preferDir(respNeuronsMB<=sign);
            RespMB{i} = max(NeuronValsMB(respNeuronsMB<=sign,:,4),[],2)';
            NeuronNMB{i} = find(respNeuronsMB<=sign);
            %RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
            pvalsMB{i} = respNeuronsMB(respNeuronsMB<=sign);
            %SpatialMB{i} = spaTuningMB(respNeuronsMB<sign);
%             OSIi = load(sprintf('Orientation-Tuning-Index-%s',NP.recordingName)).L;
%             DSIi = load(sprintf('Direction-Selectivity-Index-%s',NP.recordingName)).DSI;
%             OSI{i} = OSIi(respNeuronsMB<=sign);
%             DSI{i} = DSIi(respNeuronsMB<=sign);
            
            if noEyeMoves
                try
                    EntropMB{i} = load(sprintf('NEM-Entropies-MB-RF-respU-%s-%s.mat','0.005',NP.recordingName)).entropies;
                catch
                    EntropMB{i} = [];
                end
            else
                try
                    EntropMB{i} = load(sprintf('Entropies-MB-RF-respU-%s-%s.mat','0.005',NP.recordingName)).entropies;
                catch
                    EntropMB{i} = [];
                end

            end

            minCaxis = 0;%min(cell2mat(zscoreMB));
            maxCaxis = 1; %max(cell2mat(zscoreMB));
            %DepthUnit{i} =
            p = NP.convertPhySorting2tIc(NP.recordingDir);
            label = string(p.label');
            goodU = p.ic(:,label == 'good');

            if sum(respNeuronsMB<=sign)>0 && BrainPlots
                if DSIplot
                    colors = repmat([1 1 1],size(goodU(:,respNeuronsMB<=sign),2),1)-DSIi(respNeuronsMB<=sign); %%All black 1-DSI
                end

                if Zscoreplot
                    ncolor = 64;
                    maxZscore = 15; %by looking at z-score graph
                    possibleColors = linspace(0,maxZscore,ncolor);
                    ZscoreNeuronsMB(ZscoreNeuronsMB>maxZscore) = maxZscore;
                    for n = 1:size(goodU,2)
                        [minZS idxC(n)] = min(abs(possibleColors - ZscoreNeuronsMB(n)));
                    end
                    colors = repmat([1 1 1],size(goodU(:,respNeuronsMB<=sign),2),1)-idxC(respNeuronsMB<=sign)'/ncolor; %%All black, highest z-score
                end

                if j ==1
                    [Coor, Fig] = neuronLocation(NP,data(ex,:),goodU(:,respNeuronsMB<=sign),1,0,colors, [], 1);

                else

                    [Coor,~] = neuronLocation(NP,data(ex,:),goodU(:,respNeuronsMB<=sign),1,0,colors, Fig, 1);

                end
                j =j+1;

                unitCoorX{i} = Coor{1,1}(2,:);
                unitCoorY{i} = Coor{1,2}(2,:);
                unitCoorZ{i} = Coor{1,3}(2,:);

            elseif sum(respNeuronsMB<=sign)>0
                [Coor,~] = neuronLocation(NP,data(ex,:),goodU(:,respNeuronsMB<=sign),0,0,[], [], 1);

                unitCoorX{i} = Coor{1,1}(2,:);
                unitCoorY{i} = Coor{1,2}(2,:);
                unitCoorZ{i} = Coor{1,3}(2,:);
            else

                unitCoorX{i} = [];
                unitCoorY{i} = [];
                unitCoorZ{i} = [];
            end            

        end

        if RGsig

            zscoreMB{i} = ZscoreNeuronsMB(respNeuronsRG<=sign);
            uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
            [preferDir dirInd] = max(tuningCurve,[],2);
            PreferAngle{i} = preferDir;
            RespMB{i} = max(NeuronValsMB(respNeuronsRG<=sign,:,4),[],2)';
            NeuronNMB{i} = find(respNeuronsRG<=sign);
            %RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
            pvalsMB{i} = respNeuronsRG(respNeuronsRG<=sign);
            SpatialMB{i} = spaTuningMB(respNeuronsRG<=sign);

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

        if eachSign

            zscoreMB{i} = ZscoreNeuronsMB(respNeuronsMB<=sign);
            RespMB{i} = max(NeuronValsMB(respNeuronsMB<=sign,:,4),[],2)';
            pvalsMB{i} = respNeuronsMB(respNeuronsMB<=sign);
            NeuronNMB{i} = find(respNeuronsMB<=sign);
            nTMB = nTMB+length(ZscoreNeuronsMB(respNeuronsMB<=sign));

        end

        if bothSig

            zscoreMB{i} = ZscoreNeuronsMB(respNeuronsRG<=sign & respNeuronsMB<=sign);
            uDir = rad2deg(unique(squeeze(NeuronValsMB(:,:,5))));
            [preferDir dirInd] = max(tuningCurve,[],2);
            PreferAngle{i} = preferDir;
            RespMB{i} = max(NeuronValsMB(respNeuronsRG<=sign & respNeuronsMB<=sign,:,4),[],2)';
            NeuronNMB{i} = find(respNeuronsRG<=0.005 & respNeuronsMB<=0.005);
            %RespMB{i} = max(NeuronValsMB(:,:,4),[],2)';
            pvalsMB{i} = respNeuronsRG(respNeuronsRG<=sign & respNeuronsMB<=sign);
            SpatialMB{i} = spaTuningMB(respNeuronsRG<=sign & respNeuronsMB<=sign);


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

     %%%%%%%%%%%%%%%%%%% RG


    if contains(data.VS_ordered(ex),"RG") &&  contains(data.VS_ordered(ex),"MB")
%         RGsig =1;
%         MBsig = 0;
        NeuronValsRG= load(sprintf('RectGrid-NeuronRespCat-%s',NP.recordingName)).NeuronVals;
%        spaTuningRG = load(sprintf('Spatial-Tuning_index-RG-%s',NP.recordingName)).SpaTuning;
        respNeuronsRG = load(sprintf('RectGrid-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        ZscoreNeuronsRG = load(sprintf('RectGrid-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreU;
        
        if MBsigAnd
            
            RespRG{i} = (max(NeuronValsRG(respNeuronsMB<=sign & respNeuronsRG<=sign,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG(respNeuronsMB<=sign & respNeuronsRG<=sign);
            
%            SpatialRG{i} = spaTuningRG(respNeuronsMB<sign);
            NeuronNRG{i} = find(respNeuronsMB<=sign & respNeuronsRG<=sign);
            zscoreRG{i} = ZscoreNeuronsRG(respNeuronsMB<=sign & respNeuronsRG<=sign);

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

        end


        if MBsig
            
            RespRG{i} = (max(NeuronValsRG(respNeuronsMB<=sign,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG(respNeuronsMB<=sign & respNeuronsRG<=sign);
            
%            SpatialRG{i} = spaTuningRG(respNeuronsMB<sign);
            NeuronNRG{i} = find(respNeuronsMB<=sign );
            zscoreRG{i} = ZscoreNeuronsRG(respNeuronsMB<=sign);

            if length(ZscoreNeuronsMB(respNeuronsMB<=sign)) ~=length(ZscoreNeuronsRG(respNeuronsMB<=sign))

            end

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

        end


        if RGsig
            RespRG{i} = (max(NeuronValsRG(respNeuronsRG<=sign,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG(respNeuronsRG<=sign);
%            SpatialRG{i} = spaTuningRG(respNeuronsRG<sign);
            NeuronNRG{i} = find(respNeuronsRG<=sign);
            zscoreRG{i} = ZscoreNeuronsRG(respNeuronsRG<=sign);
        end


        if allN 
            RespRG{i} = (max(NeuronValsRG(:,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG;
  %          SpatialRG{i} = spaTuningRG;
            NeuronNRG{i} = 1:length(respNeuronsRG);
            zscoreRG{i} =0;
        end

        if eachSign
            RespRG{i} = (max(NeuronValsRG(respNeuronsRG<=sign,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG(respNeuronsRG<=sign);
            %            SpatialRG{i} = spaTuningRG(respNeuronsRG<sign);
            NeuronNRG{i} = find(respNeuronsRG<=sign);
            zscoreRG{i} = ZscoreNeuronsRG(respNeuronsRG<=sign);
            nTRG = nTRG+length(respNeuronsRG(respNeuronsRG<=sign));

        end

        if bothSig
            RespRG{i} = (max(NeuronValsRG(respNeuronsRG<=sign & respNeuronsMB<=sign,:,4),[],2))';
            pvalsRG{i} = respNeuronsRG(respNeuronsRG <=sign & respNeuronsMB<=sign);
%            SpatialRG{i} = spaTuningRG(respNeuronsRG<sign & respNeuronsMB<sign);
            NeuronNRG{i} = find(respNeuronsRG<=sign & respNeuronsMB<=sign);
            zscoreRG{i} = ZscoreNeuronsRG(respNeuronsRG<=sign & respNeuronsMB<=sign);
             EntropRG{i} = entropiesRG(respNeuronsRG<=sign & respNeuronsMB<=sign);
        end

    else
        RespRG{i} = 0;
        pvalsRG{i} = 0;
 %       SpatialRG{i} = 0;
        NeuronNRG{i} = 0;
        EntropRG{i} = 0;

    end


    %%%%%%%%%%%%%%%%%%% FFF 

    if contains(data.VS_ordered(ex),"FFF")

        NeuronValsFFF= load(sprintf('%s-FFF-respVal',NP.recordingName)).respVal;
        %        spaTuningRG = load(sprintf('Spatial-Tuning_index-RG-%s',NP.recordingName)).SpaTuning;
        respNeuronsFFF = load(sprintf('FFF-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;
        ZscoreNeuronsFFF = load(sprintf('FFF-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreU;

        if eachSign
            RespFFF{i} = NeuronValsFFF(respNeuronsFFF<=sign);
            pvalsFFF{i} = respNeuronsFFF(respNeuronsFFF<=sign);
            NeuronNFFF{i} = find(respNeuronsFFF<=sign);
            zscoreFFF{i} = ZscoreNeuronsFFF(respNeuronsFFF<=sign);
            nTFFF = nTFFF+length(respNeuronsFFF(respNeuronsFFF<=sign));

        end

        if MBsig

            RespFFF{i} = NeuronValsFFF(respNeuronsMB<=sign);
            pvalsFFF{i} = respNeuronsFFF(respNeuronsMB<=sign);
            NeuronNFFF{i} = find(respNeuronsMB<=sign);
            zscoreFFF{i} = ZscoreNeuronsFFF(respNeuronsMB<=sign);
            nTFFF = nTFFF+length(respNeuronsMB(respNeuronsMB<=sign));

        end


    else
        RespFFF{i} =0;
        pvalsFFF{i} = 0;
        NeuronNFFF{i} = 0;
        zscoreFFF{i} = 0;
        respNeuronsFFF = zeros(1,length(ZscoreNeuronsMB));

    end


     %%%%%%%%%%%%%%%%%%% SDG 

    if contains(data.VS_ordered(ex),"SDG")

        %NeuronValsFFF= load(sprintf('%s-FFF-respVal',NP.recordingName)).respVal;
        %        spaTuningRG = load(sprintf('Spatial-Tuning_index-RG-%s',NP.recordingName)).SpaTuning;
        respNeuronsSDGm = load(sprintf('SDGm-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponseM;
        ZscoreNeuronsSDGm = load(sprintf('SDGm-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreUM;
        
        respNeuronsSDGs = load(sprintf('SDGs-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponseS;
        ZscoreNeuronsSDGs = load(sprintf('SDGs-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName)).ZScoreUS;

        if eachSign
            %RespSDGm{i} = NeuronValsFFF(respNeuronsFFF<sign);
            pvalsSDGm{i} = respNeuronsSDGm(respNeuronsSDGm<=sign);
            NeuronNSDGm{i} = find(respNeuronsSDGm<=sign);
            zscoreSDGm{i} = ZscoreNeuronsSDGm(respNeuronsSDGm<=sign);
            nTSDGm = nTSDGm+length(respNeuronsSDGm(respNeuronsSDGm<=sign));


           % RespSDG{i} = NeuronValsFFF(respNeuronsFFF<sign);
            pvalsSDG{i} = respNeuronsSDGs(respNeuronsSDGs<=sign);
            NeuronNSDG{i} = find(respNeuronsSDGs<=sign);
            zscoreSDGs{i} = ZscoreNeuronsSDGs(respNeuronsSDGs<=sign);
            nTSDGs = nTSDGs+length(respNeuronsSDGs(respNeuronsSDGs<=sign));

            if length(respNeuronsSDGs) ~= length(ZscoreNeuronsMB)
                    2+2
            end

        end

        if MBsigAnd

            %Moving
            pvalsSDGm{i} = respNeuronsSDGm(respNeuronsMB<=sign & respNeuronsSDGm<=sign);
            NeuronNSDGm{i} = find(respNeuronsMB<=sign & respNeuronsSDGm<=sign);
            zscoreSDGm{i} = ZscoreNeuronsSDGm(respNeuronsMB<=sign & respNeuronsSDGm<=sign);
            nTSDGm = nTSDGm+length(respNeuronsSDGm(respNeuronsMB<=sign & respNeuronsSDGm<=sign));


            %Static
            pvalsSDG{i} = respNeuronsSDGs(respNeuronsSDGs<=sign & respNeuronsMB<=sign);
            NeuronNSDG{i} = find(respNeuronsSDGs<=sign & respNeuronsMB<=sign);
            zscoreSDGs{i} = ZscoreNeuronsSDGs(respNeuronsSDGs<=sign & respNeuronsMB<=sign);
            nTSDG = nTSDGs+length(respNeuronsSDGs(respNeuronsSDGs<=sign & respNeuronsMB<=sign));


        end

         if MBsig

            %Moving
            pvalsSDGm{i} = respNeuronsSDGm(respNeuronsMB<=sign);
            NeuronNSDGm{i} = find(respNeuronsMB<=sign);
            zscoreSDGm{i} = ZscoreNeuronsSDGm(respNeuronsMB<=sign);
            nTSDGm = nTSDGm+length(respNeuronsSDGm(respNeuronsMB<=sign));


            %Static
            pvalsSDG{i} = respNeuronsSDGs(respNeuronsMB<=sign);
            NeuronNSDG{i} = find(respNeuronsMB<=sign);
            zscoreSDGs{i} = ZscoreNeuronsSDGs(respNeuronsMB<=sign);
            nTSDG = nTSDGs+length(respNeuronsSDGs(respNeuronsMB<=sign));


        end


    else
        %RespSDGm{i} = NeuronValsFFF(respNeuronsFFF<sign);
        pvalsSDGm{i} =0;
        NeuronNSDGm{i} = 0;
        zscoreSDGm{i} = 0;

        % RespSDG{i} = NeuronValsFFF(respNeuronsFFF<sign);
        pvalsSDG{i} = 0;
        NeuronNSDG{i} = 0;
        zscoreSDG{i} = 0;
     

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

    %%%Add color vectors accordying to animals
    if string(data.Animal_ID{ex}) ~= string(data.Animal_ID{ex-1}) %wont work if you start with the first animal (noisy animal)
        animal = animal+1;
        animalName{animal} = data.Animal_ID{ex};
    end

    %%%% Counts for Venn diagram

    if Venn
    onlyMB = onlyMB+length(ZscoreNeuronsMB(respNeuronsMB<=sign & respNeuronsRG > sign & respNeuronsSDGm > sign & respNeuronsSDGs > sign));
    onlySDGm = onlySDGm+length(ZscoreNeuronsMB(respNeuronsMB>sign & respNeuronsRG > sign & respNeuronsSDGm <= sign & respNeuronsSDGs > sign));
    onlyRG = onlyRG+length(ZscoreNeuronsMB(respNeuronsMB>sign & respNeuronsRG <= sign & respNeuronsSDGm > sign & respNeuronsSDGs > sign));
    onlySDGs = onlySDGs+length(ZscoreNeuronsMB(respNeuronsMB>sign & respNeuronsRG > sign & respNeuronsSDGm > sign & respNeuronsSDGs <= sign));

    SDGmAndSDGs = SDGmAndSDGs+length(ZscoreNeuronsMB((respNeuronsSDGm <= sign & respNeuronsSDGs <= sign) & respNeuronsMB>sign & respNeuronsRG > sign ));
    SDGmAndMB = SDGmAndMB+length(ZscoreNeuronsMB((respNeuronsMB<=sign & respNeuronsSDGm <= sign) & respNeuronsRG > sign & respNeuronsSDGs > sign));
    MBandSDGs = MBandSDGs+length(ZscoreNeuronsMB((respNeuronsMB<=sign & respNeuronsSDGs <= sign) & respNeuronsRG > sign & respNeuronsSDGm > sign ));
    MBandSDGmsndSDGs = MBandSDGmsndSDGs+length(ZscoreNeuronsMB((respNeuronsMB<=sign & respNeuronsSDGm <= sign & respNeuronsSDGs <= sign) & respNeuronsRG > sign));

    MBtotal = MBtotal+length(ZscoreNeuronsMB(respNeuronsMB<=sign & respNeuronsRG > sign));
    RGtotal = RGtotal+length(ZscoreNeuronsMB(respNeuronsRG<=sign));
    SDGmtotal = SDGmtotal+length(ZscoreNeuronsMB(respNeuronsSDGm<=sign & respNeuronsRG > sign));
    SDGstotal = SDGstotal+length(ZscoreNeuronsMB(respNeuronsSDGs<=sign & respNeuronsRG > sign));
    FFFtotal = FFFtotal+length(ZscoreNeuronsMB(respNeuronsFFF<=sign));


    %%%% Z-scores accordying to venn diagram %%%%
    ZSonlyMB = [ZSonlyMB ZscoreNeuronsMB(respNeuronsMB<=sign & respNeuronsRG > sign & respNeuronsSDGm > sign & respNeuronsSDGs > sign)];
    ZSonlySDGm = [ZSonlySDGm ZscoreNeuronsSDGm(respNeuronsMB>sign & respNeuronsRG > sign & respNeuronsSDGm <= sign & respNeuronsSDGs > sign)];
    ZSonlyRG = [ZSonlyRG ZscoreNeuronsRG(respNeuronsMB>sign & respNeuronsRG <= sign & respNeuronsSDGm > sign & respNeuronsSDGs > sign)];
    ZSonlySDGs = [ZSonlySDGs ZscoreNeuronsSDGs(respNeuronsMB>sign & respNeuronsRG > sign & respNeuronsSDGm > sign & respNeuronsSDGs <= sign)];


    
        MBc = [MBc ...
            zeros(1,length(ZscoreNeuronsMB(respNeuronsMB<=sign & respNeuronsRG > sign & respNeuronsSDGm > sign & respNeuronsSDGs > sign)))+ animal];
        SDGmC = [SDGmC ...
            zeros(1,length(ZscoreNeuronsSDGm(respNeuronsMB>sign & respNeuronsRG > sign & respNeuronsSDGm <= sign & respNeuronsSDGs > sign)))+ animal];
        SDGsC = [SDGsC ...
            zeros(1,length(ZscoreNeuronsSDGs(respNeuronsMB>sign & respNeuronsRG > sign & respNeuronsSDGm > sign & respNeuronsSDGs <= sign)))+ animal];
        RGc = [RGc ...
            zeros(1,length(ZscoreNeuronsRG(respNeuronsMB>sign & respNeuronsRG <= sign & respNeuronsSDGm > sign & respNeuronsSDGs > sign))) + animal];

    end





    MB = [MB ZscoreNeuronsMB(respNeuronsMB<=sign)];
    RG = [RG ZscoreNeuronsRG(respNeuronsMB<=sign)];

    if MBsig
        MBc = [MBc ...
            zeros(1,length(ZscoreNeuronsMB(respNeuronsMB<=sign)))+ animal];
    end
    

end
%% Plot venn diagram for responsive neurons
close
[h S]= venn([MBtotal,SDGmtotal,SDGstotal],[SDGmAndMB,MBandSDGs,SDGmAndSDGs,MBandSDGmsndSDGs],'FaceColor',{"#7E2F8E","#77AC30",	"#D95319"},'EdgeColor','w');
axis off;         % Hide axis
set(gca, 'Color', 'none'); % Remove axis background
set(gcf, 'Color', 'w');  % Set figure background to white

NamesStim = {'MB','MG','SG'};
values = [onlyMB, onlySDGm, onlySDGs, SDGmAndMB,MBandSDGs,SDGmAndSDGs,MBandSDGmsndSDGs];

%Now label each zone
for i = 1:7
    if i == 1
    text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2), num2str(values(i)),"FontSize",20)
     text(S.ZoneCentroid(i,1), S.ZoneCentroid(i,2)-2, sprintf('%.2f%%', (values(i)/sum(values))*100),"FontSize",10)
    else
    text(S.ZoneCentroid(i,1)-0.3, S.ZoneCentroid(i,2)+0.7, num2str(values(i)))
     text(S.ZoneCentroid(i,1)-0.3, S.ZoneCentroid(i,2)-0.5, sprintf('%.2f%%', (values(i)/sum(values))*100),"FontSize",7)
    end
   
end

text(S.ZoneCentroid(1,1)-5, S.ZoneCentroid(1,2)+15, {NamesStim{1}, sprintf('N = %d',MBtotal)}, 'HorizontalAlignment', 'center', 'Color',"#7E2F8E")
text(S.ZoneCentroid(2,1)-2, S.ZoneCentroid(2,2)-9, {NamesStim{2}, sprintf('N = %d',SDGmtotal)},'HorizontalAlignment', 'center','Color',"#77AC30")
text(S.ZoneCentroid(3,1)-1, S.ZoneCentroid(3,2)+7, {NamesStim{3}, sprintf('N = %d',SDGstotal)},'HorizontalAlignment', 'center','Color',"#D95319")

cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
print(gcf, 'VennDiagram3groups.pdf', '-dpdf', '-r300', '-vector');

%% Plot Z-scores of neurons responding exclusively to each stimulus (accordying to venn2 diagram)

%%%%%%%%% Prepare data for summary plots

%%%%Add SDGm, and SDGs (4 color vectors) 

colormapUsed = parula(animal).*0.8;

% Define x-values for each group
xMB = ones(size(ZSonlyMB))' * 1;     % Group 1 at x = 1
xSDGm = ones(size(ZSonlySDGm))' * 2; % Group 2 at x = 2
xSDGs = ones(size(ZSonlySDGs))' * 3; % Group 3 at x = 3
xRG = ones(size(ZSonlyRG))' * 4;     % Group 4 at x = 4

% Combine data
x = [xMB; xSDGm; xSDGs; xRG];
y = log([ZSonlyMB';ZSonlySDGm';ZSonlySDGs';ZSonlyRG']);


% Combine all color indices
allColorIndices = [MBc, SDGmC, SDGsC, RGc];

% Calculate the number of points for each color (this is the actual count)
colorCounts = histcounts(allColorIndices, 1:n+1); % Count occurrences of each color

f = figure;
% ---- Swarmchart (Larger Left Subplot) ----
subplot(1, 10, 1:7); % Takes most of the space
swarmchart(x, y, 20, colormapUsed(mod(allColorIndices-1, n) + 1, :), 'filled'); % Marker size 50
xticks([1 2 3 4]);
xticklabels({'MB', 'SDGm', 'SDGs', 'SB'});
xlabel('Stimuli');
ylabel('ln(Z-score)');
grid on;
set(f,'Color','w')

% ---- Narrower Vertical Rectangle (Smaller Right Subplot) ----
% Change subplot size and position for the rectangle
subplot(1, 10, 8:9); % Much narrower and closer to the swarm chart
hold on;

totalPoints = sum(colorCounts); % Total points across all colors
yBottom = 0;
rectWidth = 0.3; % Even narrower for better fit

% Plot each section of the rectangle and add color index text
for colorIdx = 1:animal
    count = colorCounts(colorIdx); % Get the count for this color
    height = count / totalPoints;  % Proportional height of each color category
    
    % Plot the rectangle for this color category
    rectangle('Position', [0, yBottom, rectWidth, height], 'FaceColor', colormapUsed(colorIdx, :), 'EdgeColor', 'k');
    
    % Add color index (vertical text) next to each section
    text(rectWidth + 0.05, yBottom + height/2, animalName{colorIdx}, ...
        'Rotation', 90, 'VerticalAlignment', 'middle', 'FontSize', 12);
    
    % Add point count (horizontal text) inside the section
    text(rectWidth/2, yBottom + height/2, num2str(count), ...
        'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color','w');
    
    yBottom = yBottom + height; % Move up for the next section
end

% Adjust axes
axis off;
xlim([0, 0.4]); % Keeps it narrow but closer
ylim([0, 1]);

title('N. per Animal');
hold off;


cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
print(gcf, 'Z-score-MB-SDG-RG_Exclusive neurons.pdf', '-dpdf', '-r300', '-vector');
%% Plot Z-scores of neurons responding to moving ball 






%% plot colorbar for brain plot
c= colorbar;
colormap(repmat([1 1 1],64,1)-repmat(linspace(0,1,64),3,1)')
c.Title.String = "DSI";

print('DSI_index_brain','-dpdf', '-r1000');

%% plot variables against X,y and Z unit position
[uniqueStrings,~,indexAnimal] = unique([animalID{:}],'stable');
color1 = [];
for i = 1:length(animalID) %each color is an animal.
    
    color1 = [color1 zeros(1,length(RespMB{i}))+indexAnimal(i)];


end

colormap('jet'); % Use a colormap (e.g., 'jet', 'parula', etc.)
colors = colormap; % Get colormap colors
n_colors = size(colors, 1); % Number of colors in the colormap

% Normalize indices to match the colormap range
mapped_colors = colors(ceil((color1 / max(color1)) * n_colors), :);


DSIm = cell2mat(DSI')';
entropM = cell2mat(EntropMB);
X = cell2mat(unitCoorX);
Y = cell2mat(unitCoorY);
Z = cell2mat(unitCoorZ);


y = DSIm;
fig = tiledlayout(3,3);

nexttile
x = X-min(X);
scatter(x,y,50,mapped_colors,'filled')
xlabel('Medial to Lateral')
ylabel('DSI')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));

nexttile
x = Y-min(Y);
scatter(x,y,50,mapped_colors,'filled')
xlabel('Rostral to Caudal')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));

nexttile
x = Z;
scatter(x,y,50,mapped_colors,'filled')
xlabel('Depth')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));


y = cell2mat(zscoreMB);
nexttile
x = X-min(X);
scatter(x,y,50,mapped_colors,'filled')
xlabel('Medial to Lateral')
ylabel('Z-score')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));

nexttile
x = Y-min(Y);
scatter(x,y,50,mapped_colors,'filled')
xlabel('Rostral to Caudal')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));

nexttile
x = Z;
scatter(x,y,50,mapped_colors,'filled')
xlabel('Depth')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));

y = entropM;

nexttile
x = X-min(X);
scatter(x,y,50,mapped_colors,'filled')
xlabel('Medial to Lateral')
ylabel('Entropy')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));

nexttile
x = Y-min(Y);
scatter(x,y,50,mapped_colors,'filled')
xlabel('Rostral to Caudal')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));

nexttile
x = Z;
scatter(x,y,50,mapped_colors,'filled')
xlabel('Depth')
mdl = fitlm(x, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));




%% Plot regresion line for values of neurons across x,y,and z 
figure;
x = DSIm;
y = cell2mat(zscoreMB);
scatter(x,y,50,mapped_colors,'filled')
xlabel('DSI')
ylabel('Z-score')
mdl = fitlm(DSIm, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));


figure;
x = DSIm;
y = entropM;
scatter(x,y,50,mapped_colors,'filled')
xlabel('DSI')
ylabel('Entropy')
mdl = fitlm(DSIm, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));



figure;
x = cell2mat(zscoreMB);
y = entropM;
scatter(x,y,50,mapped_colors,'filled')
xlabel('Z-score')
ylabel('Entropy')
mdl = fitlm(DSIm, y);
x_fit = linspace(min(x), max(x), 100); % Generate x values for the line
y_fit = predict(mdl, x_fit');          % Predict corresponding y values
hold on;plot(x_fit, y_fit, 'r-', 'LineWidth', 2); % Red fitted line
R_squared = mdl.Rsquared.Ordinary; % Extract R^2 value
title(sprintf('Rsquared = %d',R_squared));

figure; scatter(1:max(indexAnimal),ones(1,max(indexAnimal)),50,colors(ceil((unique(indexAnimal) / max(indexAnimal)) * n_colors), :))

%% Prepare data for summary plots
MB = cell2mat(RespMB);
RG = cell2mat(RespRG);

%%%%Add SDGm, and SDGs

[uniqueStrings,~,indexAnimal] = unique([animalID{:}], 'stable');
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





%% Plot Z-scores (for neurons responsive to MB) MBcompa
MB = cell2mat(zscoreMB)';
RG = cell2mat(zscoreRG)';
SDGm = cell2mat(zscoreSDGm)';
SDGs = cell2mat(zscoreSDGs)';

colormapUsed = parula(animal).*0.8;

% eMB = cell2mat(EntropMB);
% eRG = cell2mat(EntropRG);

%y =log10([MB; RG; SDGm; SDGs; MB-RG; MB-SDGm; MB-SDGs]);

y =([MB; RG; MB-RG]);


x =[];
% Combine data
for i =1:length(y)/length(MB)
    x = [x;ones(size(MB))*i];
end

% Combine all color indices
allColorIndices = repmat(MBc',length(y)/length(MB),1);

% Calculate the number of points for each color (this is the actual count)
colorCounts = histcounts(MBc', 1:n+1); % Count occurrences of each color

f = figure;
% ---- Swarmchart (Larger Left Subplot) ----
subplot(1, 15, 1:13); % Takes most of the space
swarmchart(x, y, 5, colormapUsed(mod(allColorIndices-1, n) + 1, :), 'filled','MarkerFaceAlpha',0.7); % Marker size 50
xticks(1:7);
xticklabels({'MB', 'SB', 'MB-SB'});
ylabel('Z-score');
set(f,'Color','w')
%set(gca, 'YScale', 'log')
yline([0],'LineWidth',2)
ylim([-5 20])
f.Position = [680   680   359   298];

% % ---- Narrower Vertical Rectangle (Smaller Right Subplot) ---- to plot
% neuron percentages
% % Change subplot size and position for the rectangle
% subplot(1, 15, 14:15); % Much narrower and closer to the swarm chart
% hold on;
% 
% totalPoints = sum(colorCounts); % Total points across all colors
% yBottom = 0;
% rectWidth = 0.3; % Even narrower for better fit
% 
% % Plot each section of the rectangle and add color index text
% for colorIdx = 1:animal
%     count = colorCounts(colorIdx); % Get the count for this color
%     height = count / totalPoints;  % Proportional height of each color category
%     
%     % Plot the rectangle for this color category
%     rectangle('Position', [0, yBottom, rectWidth, height], 'FaceColor', colormapUsed(colorIdx, :), 'EdgeColor', 'k');
%     
%     % Add color index (vertical text) next to each section
%     text(rectWidth + 0.05, yBottom + height/2, animalName{colorIdx}, ...
%         'Rotation', 90, 'VerticalAlignment', 'middle', 'FontSize', 12);
%     
%     % Add point count (horizontal text) inside the section
%     text(rectWidth/2, yBottom + height/2, num2str(count), ...
%         'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', 'FontSize', 12, 'Color','w');
%     
%     yBottom = yBottom + height; % Move up for the next section
% end
% 
% % Adjust axes
% axis off;
% xlim([0, 0.4]); % Keeps it narrow but closer
% ylim([0, 1]);
% 
% title('N. per Animal');
% hold off;


cd('G:\My Drive\2024_Large scale mapping of visual responses in the ventral pallium reveal tuning to visual motion\Figure\Fig1_panels')
print(gcf, 'Z-score-MB-RG_MovBall_resp_neurons.pdf', '-dpdf', '-r300', '-vector');



%%%Old method:
% %values = ([eMB eRG eMB-eRG]);
% 
% % %%randomly select 10 neurons but from animals that have full receptive
% % %%field tuning == PV35 & PV139 19:20 40:43
% % %goodRFr = cell2mat(zscoreMB(1))-cell2mat(zscoreRG(19:24));
% % goodRFr = cell2mat(EntropMB(1))-cell2mat(EntropRG(1));
% % %MBh = length(MB)+length(RG)+length(cell2mat(zscoreMB(1)))+find(goodRFr<0);
% % MBh = length(MB)+length(RG)+find(goodRFr<0);
% % RGh = length(MB)+length(RG)+find(goodRFr>0);
% % 
% % seed1 = 42;
% % 
% % seed2 = 43;
% % randN = 5;
% % rng(42);
% % rMB = 1 + round((length(MBh) - 1) * rand(1, randN));
% % rng(42);
% % rRG = 1 + round((length(RGh) - 1) * rand(1, randN));
% % 
% % %selected neurons:
% 
% %selecN = {[tablePosN((MBh(rMB)-length(MB)*2));NeuronID((MBh(rMB)-length(MB)*2))];...
%   %  [tablePosN(RGh(rRG)-length(MB));NeuronID(RGh(rRG)-length(MB))]};
% 
% 
% %randomNeur = 1 + round((length(cats) - 1) * rand(1, 10));
% 
% T = table(cats',values',color','VariableNames',{'cats','values','color'});
% figure;
% % 
% s = swarmchart(T.cats, T.values, 10, T.color,'filled','MarkerFaceAlpha',0.5);
% % 
% % hold on %%Chose before, 2 and 5;
% % chooseNb = 1;
% % swarmchart(T.cats([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]), ...
% %     T.values([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]), 70,...
% %     T.color([MBh(rMB(chooseNb))-length(MB)*2 MBh(rMB(chooseNb))-length(MB) MBh(rMB(chooseNb))]),...
% %     'filled','MarkerEdgeColor','r','LineWidth',2);
% 
% % hold on
% % chooseNr = 1;
% % swarmchart(T.cats([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]),...
% %     T.values([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]), 50, ...
% %     T.color([RGh(rRG(chooseNr))-length(MB)*2 RGh(rRG(chooseNr))-length(MB) RGh(rRG(chooseNr))]),...
% %     'filled','MarkerEdgeColor','g','LineWidth',2);
% % 
% % 
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
% 
% 
% % myColormap = [
% %     1, 0, 0;    % Red
% %     0, 1, 0;    % Green
% %     0, 0, 1;    % Blue
% %     1, 0, 1;    % Magenta
% %     0, 0, 0;    % Black
% % ]*0.7;
% % colormap(gca,myColormap)
% 
% % myColormap = [
% %     0, 0, 1;    % Blue
% %     1, 0, 1;    % Magenta
% %     0, 0, 0;    % Black
% % ]*0.7;
% % colormap(gca,myColormap)
% 
% 
% 
% % hold on;
% % uniqueGroups = unique(T.color);  % Get unique group values
% %  % Get the same colormap used in the plot
% % 
% % for i = [1 3]
% %     % Plot invisible points for legend creation
% %     scatter(nan, nan, 100, colors(i,:), 'filled');
% %     
% % end
% % 
% % 
% % 
% hold on
% for i =1:3
% 
%     groupMean = mean(T.values(T.cats==categorical(i))); 
% 
%      plot([i-0.5, i+0.5], [groupMean, groupMean], 'r-', 'LineWidth', 2);  % Short horizontal line for the mean
% 
% end
% % 
% % legend({'PV35', 'PV103', 'SA5'}, 'Location', 'best');
% grid on
% set(gcf,'Color','w');%
% xticklabels({'Moving','Static','Mov. - Stat.'})
% ylabel('Z-score')
% yL = ylim;
% %ylim([yL(1) 60])
% cd('\\sil3\data\Large_scale_mapping_NP\lizards\SummaryFigs')
% %c = colorbar;
% % title(c, 'Animals (N=5)');%caxis([1 5])
% %c.Ticks = [linspace(1.5,4.5,5)];
% % c.Ticks = [linspace(1.5,2.5,animalN)];
% % c.TickLabels = {'PV67','PV139',  'PV35'};
% 
% %print(gcf,'NEM-AWAKESIg(MB)-entropy-MB-RG(red-High-MB)-(green-High-RG)-Color-Animal_1N)','-dpng')
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

