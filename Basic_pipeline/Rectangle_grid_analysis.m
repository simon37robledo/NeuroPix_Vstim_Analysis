%%% Rect Grid Analysis


cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%Optionall
plotRasters =1;
heatMap = 0;
ex=1;


%%
for ex =44%1:size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
    path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    cd(path)
    NP = NPAPRecording(path);


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


    if isempty(rectFiles)
        %disp()
        w= sprintf('No rectangle grid files where found in %s. Skipping into next experiment.',NP.recordingName);
        warning(w)
        continue
    end

    j =1;
    
    %%%Get only rectangle grid files, no novelty. 
    %Extract numbers (date) and sort the cell array
    numbers = cellfun(@(str) str2double(regexp(str, '(?:[^_]*_){4}(\d+)_(\d+)', 'tokens', 'once')), rectFiles, 'UniformOutput', false);


    % Convert to a matrix for sorting
    numbers = cell2mat(numbers);
    j=1;
    combinedNumbers = zeros(1,numel(rectFiles));
    for n = 1:2:numel(rectFiles)*2
        combinedNumbers(j) = numbers(:, n) * 1000 + numbers(:, n+1);
        j=j+1;

    end

    % Sort based on the 4th and then 5th number
    [~, sortIdx] = sort(combinedNumbers); % Sort by 4th, then 5th number
    rectFiles = rectFiles(sortIdx); % Sort the cell array by time in file

    VSordered = strsplit(data.VS_ordered{ex},',');
    RGpos = find(VSordered=="RG");
    OBpos = find(VSordered=="OB");
    OBCpos = find(VSordered=="OBC");

    [orderVS orderVSIndex] = sort([RGpos OBpos OBCpos]);

    selecFiles = rectFiles(orderVSIndex(1));


    if size(rectFiles) ~= [0 0]

        for i = selecFiles
            rect= load(stimDir+"\"+string(i));

           
                %%New exp
                seqMatrix = [seqMatrix cell2mat(rect.VSMetaData.allPropVal(find(strcmp(rect.VSMetaData.allPropName,'pos'))))];
                if ex >18 || ex <28
                    sizes = [sizes cell2mat(rect.VSMetaData.allPropVal(find(strcmp(rect.VSMetaData.allPropName,'tilingRatios'))))];
                end
                stimDurStats = cell2mat(rect.VSMetaData.allPropVal(42))*1000;
                interStimStats = cell2mat(rect.VSMetaData.allPropVal(32))*1000;
     
            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    positionsMatrix = [cell2mat(rect.VSMetaData.allPropVal(find(strcmp(rect.VSMetaData.allPropName,'pos2X'))))...
            ,cell2mat(rect.VSMetaData.allPropVal(find(strcmp(rect.VSMetaData.allPropName,'pos2Y'))))];%NewExp
    rectSizes = cell2mat(rect.VSMetaData.allPropVal(find(strcmp(rect.VSMetaData.allPropName,'rectSide'))));

    if ex<19 || ex >27
        %OldExp
        sizes = repmat(cell2mat(rect.VSMetaData.allPropVal(5)),1,length(seqMatrix));
    end

    nSize = length(unique(sizes));
    % Triggers:

    Ordered_stims= strsplit(data.VS_ordered{ex},',');
    containsRG = cellfun(@(x) contains(x,'RG'),Ordered_stims);
    ttlInd = find(containsRG);

    [stimOn stimOff] = NPdiodeExtract(NP,newDiode,0,"RG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff] = NPdiodeExtract(NP,0,0,"RG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A


    %missedT = find((stimOff-stimOn)<500); %Missed tri als

    stimDur = mean(stimOff'-stimOn');
    interStim = mean(stimOn(2:end)-stimOff(1:end-1));

    A = [stimOn seqMatrix' sizes'];

    [C indRG]= sortrows(A,[2 3]);

    %Sort directions:

    directimesSorted = C(:,1)';

    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

    %Good units

    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

    p = NP.convertPhySorting2tIc(NP.recordingDir);

    %Select good units
    label = string(p.label');

    goodU = p.ic(:,label == 'good');

    % Build raster plots per unit

    if plotRasters == 1
        for plotRasters=1
            cd(NP.recordingDir)
            if ~exist(path+"\Figs",'dir')
                mkdir Figs
            end
            cd(NP.recordingDir + "\Figs")

            bin=20;
            win=stimDur+interStim;
            preBase = round(interStim/20)*10;

            [Mr]=BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(win/bin));

            [nT,nN,nB] = size(Mr);
            %indRG --> sorted infexes

            trialsPerCath = length(seqMatrix)/(length(unique(seqMatrix)));

            PositionsTotal = positionsMatrix(seqMatrix(indRG),:);


            [posS,indexX] = sortrows(PositionsTotal,1); %Sort first dimension because tile layout moves through columns

            MrS = Mr(indexX,:,:);



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

        end
    end

    if heatMap ==1
    % Build heat map

    cd(NP.recordingDir)
    if ~exist(path+"\Figs",'dir')
        mkdir Figs
    end
    cd(NP.recordingDir + "\Figs")

    trialDiv  = length(seqMatrix)/length(unique(seqMatrix))/nSize;

    offsetR=400; %ms for offset response

    bin =1;
    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round(stimDur)/bin);
    [Mro] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted+stimDur)/bin),round(interStim)/bin);

    [nT,nN,NB] = size(Mr);
    [nTo,nNo,NBo] = size(Mro);


    MrC = zeros(round(nT/trialDiv),nN, NB+NBo);

    %%Create summary of identical trials

    for u = 1:length(goodU)
        j=1;

        for i = 1:trialDiv:nT

            meanRon = mean(squeeze(Mr(i:i+trialDiv-1,u,:)));

            meanRoff =  mean(squeeze(Mro(i:i+trialDiv-1,u,:)));

            MrC(j,u,:) = [meanRon meanRoff]; %Combine on and off response

            j = j+1;

        end
    end
   

    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin));

    Nb2=  mean(Mb,3);

    Nbase = mean(Mb,[1 3]);

    MrMean = mean(MrC,3)-Nbase;

    cd(NP.recordingDir + "\Figs")
%     respU = [138,134,127,124,123,121,118,112]; %PV67-1 %check Offset response also.
%     respU = [51,28];
    % Create a meshgrid of coordinates

    reduF = 10; 
    screenSide = rect.VSMetaData.allPropVal{find(strcmp(rect.VSMetaData.allPropName,'rect'))}; %Same as moving ball

    screenRed = screenSide(4)/reduF;
    [x, y] = meshgrid(1:screenRed, 1:screenRed);
    %screenSide = (max(rect.VSMetaData.allPropVal{21,1}.Y4{1,4},[],'all'))

    RFuT = zeros(screenRed,screenRed,length(u));

    j=1;

    %u =12;

    %  RFu = zeros(screenSide,screenSide,length(C)/);

    pxyScreen = zeros(screenRed,screenRed);

    VideoScreen = zeros(screenRed,screenRed,length(C)/trialDiv);

    prop = find(strcmp(rect.VSMetaData.allPropName,'rectData'));


    for i = 1:trialDiv:length(C)

        xyScreen = zeros(screenRed,screenRed)'; %%Make calculations if sizes>1 and if experiment is new and the shape is a circle.

        if nSize>1 && string(rect.VSMetaData.allPropVal{7}) == "circle"
            
            Xc = round((rect.VSMetaData.allPropVal{prop,1}.X2{1,C(i,3)}(C(i,2))-rect.VSMetaData.allPropVal{prop,1}.X1{1,C(i,3)}(C(i,2)))/2)+rect.VSMetaData.allPropVal{prop,1}.X1{1,C(i,3)}(C(i,2));%...
            %-min(rect.VSMetaData.allPropVal{21,1}.Y2{1,4},[],'all'))+conversion;
            Xc = Xc/reduF;
            
            Yc = round((rect.VSMetaData.allPropVal{prop,1}.Y4{1,C(i,3)}(C(i,2))-rect.VSMetaData.allPropVal{prop,1}.Y1{1,C(i,3)}(C(i,2)))/2)+rect.VSMetaData.allPropVal{prop,1}.Y1{1,C(i,3)}(C(i,2));%...
            %-min(rect.VSMetaData.allPropVal{21,1}.Y2{1,4},[],'all'))+conversion;
            Yc = Yc/reduF;
           
            r = round((rect.VSMetaData.allPropVal{prop,1}.X2{1,C(i,3)}(C(i,2))-rect.VSMetaData.allPropVal{prop,1}.X1{1,C(i,3)}(C(i,2)))/2);
            r= r/reduF;
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

    VD = repmat(VideoScreen,[1 1 1 nN]);

    Res = reshape(MrMean,[1,1,size(MrMean,1),nN]);

    RFu = squeeze(sum(VD.*Res,3));

    figure;imagesc(normRFu(:,:,51))

    normMatrix = repmat(pxyScreen,[1,1,nN]).*reshape(Nbase,[1,1,nN]);

    normRFu = RFu./normMatrix;

    save(sprintf('RFuStatic-%s',NP.recordingName),"normRFu")





%     epsilon = 0.01;
% 
%     denom = mad(Nb2(u),0)+epsilon; %mean(Mb,0)+epsilon; %
% 
%     % mSpk = mean(spkRateR);
% 
% %     M = ( MrMean(:,u)' - (Nbase(u) + MrMean(:,u)')./2)./denom;
% 
%     RFu = mean(bsxfun(@times, VideoScreen, reshape(M, 1, 1, [])),3);
%     %figure;imagesc(RFu)
%     fig = figure;imagesc(RFu./pxyScreen)
%     caxis([0 0.02]);
%     colorbar; max(RFu./pxyScreen,[],'all')
%     xlabel('X pixels')
%     ylabel('Y pixels')
%     title(sprintf('RFu-%d-Static',u))
%     prettify_plot
% 
%     print(fig,sprintf('%s-RFu-%d-Static.png',NP.recordingName,u),'-dpng')
%     close

    end
end












