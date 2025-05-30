%%% Rect Grid Analysis


cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%Optionall
bombcelled = 0;
plotRasters =0;
heatMap = 1;
plotHeatMap =0;
calculateEntropy =0;
ex=1;
responseDo = 1;
redoResp =0;
noEyeMoves = 1;
newDiode =0;

Shuffling_baseline=0;
repeatShuff =0;
trialThres =0.6;

%GoodRecordingsPV =[15:20,40:43];
SpatialTuning = 0;

Awake = [44:48];


%GoodRecordingsPV =[1:21,40:43,49:54];

%%
for ex = 51 %selecN{1}(1,:)%20%GoodRecordingsPV%[1:20,28:32,40:48]%1:size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
     %1. Load NP class
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


    if ~contains(data.VS_ordered(ex),"RG")
        %disp()
        w= sprintf('No rectangle grid files where found in %s. Skipping into next experiment.',NP.recordingName);
        warning(w)
        continue
    end

    
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

    selecFiles = rectFiles(orderVSIndex(1)); %Select Rectangle grid files


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
    nPos = length(unique(seqMatrix));

    trialDivision = length(seqMatrix)/nSize/nPos;
    % Triggers:

    Ordered_stims= strsplit(data.VS_ordered{ex},',');
    containsRG = cellfun(@(x) contains(x,'RG'),Ordered_stims);
    ttlInd = find(containsRG);

    [stimOn stimOff] = NPdiodeExtract(NP,newDiode,0,"RG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff] = NPdiodeExtract(NP,0,0,"RG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A


    %missedT = find((stimOff-stimOn)<500); %Missed trials

    stimDur = mean(stimOff'-stimOn');
    stimInter = mean(stimOn(2:end)-stimOff(1:end-1));

    preBase = round(stimInter/2);

    A = [stimOn seqMatrix' sizes'];

    [C indRG]= sortrows(A,[2 3]);

    %Sort directions:

    directimesSorted = C(:,1)';

    stims = stimOn';

    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

    %Good units

    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

    if bombcelled
        p = NP.convertPhySorting2tIc(NP.recordingDir,0,1);
    else
        p = NP.convertPhySorting2tIc(NP.recordingDir);
    end

    %Select good units
    label = string(p.label');

    goodU = p.ic(:,label == 'good');

    if isempty(goodU)
        %disp()
        w= sprintf('No somatic neurons in %s. Skipping into next experiment.',NP.recordingName);
        warning(w)
        continue
    end

    bin =1;

    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

    Mb = mean(Mb,3); %mean across time bins

    spkRateBM = mean(Mb);

    if responseDo ==1
        %4. Create window to scan rasters and get the maximum response
        duration =300; %window in ms, same as in MB
        Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((stimDur+duration)/bin)); %response matrix
        [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');

        [nT,nN,nB] = size(MrC);


        window_size = [1, round(duration/bin)];

        %5. Initialize the maximum mean value and its position


        max_position = zeros(nN,2);
        max_mean_value = zeros(1,nN);
        max_mean_valueB = zeros(1,nN);


        NeuronVals = zeros(nN,nT/trialDivision,7); %Each neuron, which has a matrix where the first column is maxVal of bins, 2nd and 3rd position of window in matrix...
        % 4th Z-score.
        % responsive compared to the baseline.

        %[1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 17, 20, 22, 23, 24, 28, 32, 34, 36]

        if ~isfile(sprintf(sprintf('RectGrid-NeuronRespCat-%s.mat',NP.recordingName))) || redoResp

            mergeTrials = trialDivision;

            %Baseline = size window

            [Mbd] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

            %Merge trials:

           
            Bd = reshape(Mbd, [mergeTrials, size(Mbd, 1)/mergeTrials, size(Mbd, 2), size(Mbd,3)]);

            Mbd2 = squeeze(mean(Bd, 1));

            %Merge trials:

            mergeTrials = trialDivision;

            B = reshape(MrC, [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);

            % Take the mean across the first dimension (rows)
            Mr2 = squeeze(mean(B, 1));

            [nT,nN,nB] = size(Mr2);
            cd(NP.recordingDir)
            
            

            %Real data:
            %A = [stimOn seqMatrix' sizes'];
            %
            %[C indRG]= sortrows(A,[2 3]);

            for u =1:nN
                % Slide the window over the matrix
                %unit matrix
                max_mean_value(u) = -Inf; %General max? not needed
                max_mean_valueB(u)=-Inf;
                NeuronRespProfile = zeros(nT,7); %4 columns plus: ofsett, dir, size, speed, frec.

                k =1;
                max_position_Trial = zeros(nT,2); %Create 2 matrices, for mean inside max window, and for position of window, for each trial category
                max_mean_value_Trial = zeros(1,nT);
                max_mean_value_TrialB = zeros(1,nT);

                for i = 1:nT %Iterate across trials
                    uM = squeeze(Mr2(i,u,:))';%Create matrix per unique trial conditions
                    uMB = squeeze(Mbd2(i,u,:))';%Create matrix per unique trial conditions

                    max_mean_value_Trial(k) = -Inf;
                    max_mean_value_TrialB(k) = -Inf;

                    for j = 1:2:size(uM, 2) - window_size(2) + 1 %Iterate across bins
                        % Extract the sub-matrix
                        sub_matrix = uM(j:min(j+window_size(2)-1,end)); %Create matrix of size window per bin
                        sub_matrixB = uMB(j:min(j+window_size(2)-1,end));

                        % Compute the mean value
                        mean_value = mean(sub_matrix(:)); %Compute mean of each window
                        mean_valueB = mean(sub_matrixB(:)); %Compute mean of each window

                        % Update the maximum mean value and its position (a
                        % window is selected across each trial division
                        if mean_value >  max_mean_value_Trial(k)
                            max_mean_value_Trial(k) = mean_value;
                            max_position_Trial(k,:) = [i j];
                        end

                        if mean_valueB >  max_mean_value_TrialB(k)
                            max_mean_value_TrialB(k) = mean_valueB;
                        end

                    end
                    %Save across each trial (in a row) the max bin window
                    %1%. Response
                    NeuronRespProfile(k,1) = max_mean_value_Trial(k);

                    %2%. WindowTrial
                    NeuronRespProfile(k,2) = max_position_Trial(k,1);

                    %3%. WindowBin
                    NeuronRespProfile(k,3) = max_position_Trial(k,2);

                    %4%. Resp - Baseline
                    % NeuronRespProfile(i,4) = (max_mean_value_Trial(i) - (spkRateBM(u)+max_mean_value_Trial(i))/2)/denom(u); %Zscore
                    %NeuronRespProfile(k,4) = (max_mean_value_Trial(k) - spkRateBM(u))/denom(u); %Zscore
                    NeuronRespProfile(k,4) = max_mean_value_Trial(k)-max_mean_value_TrialB(k);
                    maxWindow = squeeze(Mr(k*trialDivision-trialDivision+1:k*trialDivision,u,max_position_Trial(k,2):max_position_Trial(k,2)+duration-1));
                    emptyRows = sum(all(maxWindow == 0, 2));
                    if emptyRows/trialDivision > 0.6
                        NeuronRespProfile(k,7) = 0;
                    else
                        NeuronRespProfile(k,7) = 1; %%enough trials with spiking
                    end
                    %Assign visual stats to last columns of NeuronRespProfile. Select
                    %according  to trial (d) the appropiate parameter > directions'offsets' sizes' speeds' freq'
                    NeuronRespProfile(k,5) = C(i*mergeTrials,2);
                    NeuronRespProfile(k,6) = C(i*mergeTrials,3);

                    k = k+1;

                end

                %figure;imagesc(uM);xline(max_position_Trial(i,2));xline(max_position_Trial(i,2)+window_size(2))

                NeuronVals(u,:,:) = NeuronRespProfile;
            end
            save(sprintf('RectGrid-NeuronRespCat-%s',NP.recordingName),"NeuronVals")
        else
            NeuronVals = load(sprintf('RectGrid-NeuronRespCat-%s',NP.recordingName)).NeuronVals;
        end


        %%% %Spatial tuning:

        if SpatialTuning ==1

        boot_means = load(sprintf('RectGrid-Base-Boot-1000-%s',NP.recordingName)).boot_means; 
        ZScoreU = zeros(length(unique(seqMatrix)),size(goodU,2));

        N_bootstrap = 1000;
 
        for u = 1:size(goodU,2)

            response = mean(reshape(squeeze(NeuronVals(u,:,1)),[nSize size(NeuronVals,2)/nSize])); %Take the mean across all sizes: 

            ZScoreU(:,u) = (response-mean(boot_means(:,u)))/(std(boot_means(:,u))+1/(N_bootstrap*trialDivision));

        end

        [PreferPos ind]= max(abs(ZScoreU),[],1);

        SpaTuning = zeros(1,size(goodU,2));

        for u =1:size(goodU,2)

            SpaTuning(u) = 1- mean(abs(ZScoreU(setdiff(1:size(ZScoreU,1), ind(u)),u)))./abs(PreferPos(u));
        end

        save(sprintf('Spatial-Tuning_index-RG-%s',NP.recordingName),"SpaTuning")

        end

        if Shuffling_baseline

            if ~isfile(sprintf('RectGrid-pvalsBaselineBoot-%d-%s.mat',N_bootstrap,NP.recordingName))||repeatShuff==1

                baseline = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((preBase)/bin));
                baseline = single(baseline);
                [nT,nN,nB] = size(baseline);

                % Bootstrapping settings
                N_bootstrap = 1000; % Number of bootstrap iterations
                boot_means = zeros(N_bootstrap, nN,'single');
                resampled_indicesTr = single(randi(nT, [nT, N_bootstrap]));% To store bootstrapped means
                resampled_indicesTi = single(randi(nB, [nB, N_bootstrap]));

                kernel = ones(trialDivision, duration) / (trialDivision * duration); % Normalize for mean
                % Start a parallel pool (if not already started)
                if isempty(gcp('nocreate'))
                    parpool; % Start a pool with the default number of workers
                end

                tic
                parfor i = 1:N_bootstrap
                    % Resample trials with replacement
                    resampled_trials = baseline(resampled_indicesTr(:, i), :,resampled_indicesTi(:, i));
                    for ui = 1:nN
                        % Extract the slice for the current unit (t x b matrix)
                        slice = resampled_trials(:, ui, :);
                        slice = squeeze(slice); % Result is t x b

                        % Compute the mean using 2D convolution
                        means = conv2(slice, kernel, 'valid'); % 'valid' ensures the window fits within bounds

                        % Find the maximum mean in this slice
                        boot_means(i, ui) = max(means(:));
                    end
                end
                toc

                [respVal,respVali]= max(NeuronVals(:,:,1),[],2);

                Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((stimDur+duration)/bin)); %response matrix

                %%% Calculate p-value & Filter out neurons in which max response window is empty for more than
                %%% 60% of trials

                pvalsResponse = zeros(1,nN);
                ZScoreU = zeros(1,nN);

                for u = 1:nN
                    posTr = NeuronVals(u,respVali(u),2);
                    posBin = NeuronVals(u,respVali(u),3);

                    maxWindow = squeeze(Mr(posTr*trialDivision-trialDivision+1:posTr*trialDivision,u,posBin:posBin+duration-1));

                    emptyRows = sum(all(maxWindow == 0, 2));

                    pvalsResponse(u) = mean(boot_means(:,u)>respVal(u));
                    ZScoreU(u) = (respVal(u)-mean(boot_means(:,u)))/(std(boot_means(:,u))+1/(N_bootstrap*trialDivision));

                    if emptyRows/trialDivision >= trialThres
                        pvalsResponse(u) = 1;
                    end

                end
                cd(NP.recordingDir)
                save(sprintf('RectGrid-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName),'pvalsResponse')
                save(sprintf('RectGrid-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName),'ZScoreU')
                save(sprintf('RectGrid-Base-Boot-%d-%s',N_bootstrap,NP.recordingName),'boot_means')
            end

        end


    end

    % Build raster plots per unit

    if plotRasters == 1
        for plotRasters=1
            cd(NP.recordingDir)
            if ~exist(path+"\Figs",'dir')
                mkdir Figs
            end
            cd(NP.recordingDir + "\Figs")
            
            bin=30;
            win=stimDur+stimInter;
            preBase = round(stimInter/20)*10;

            [Mr]=BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(win/bin));

            Mr = ConvBurstMatrix(Mr,fspecial('gaussian',[1 10],3),'same');

            [nT,nN,nB] = size(Mr);
            %indRG --> sorted infexes

            trialsPerCath = length(seqMatrix)/(length(unique(seqMatrix)));

            posX = squeeze(NeuronVals(:,:,3));
            posY = squeeze(NeuronVals(:,:,2));

            eNeuron = 1:length(goodU);
            eNeuron =[60];%selecN{1}(2,selecN{1}(1,:)==ex);

            [respVal respVali]= max(NeuronVals(:,:,1),[],2);


            for u =eNeuron%1:length(goodU)
                t = tiledlayout(sqrt(max(seqMatrix)), sqrt(max(seqMatrix)),'TileSpacing','compact');


             

                if trialsPerCath>20
                    j=1;
                    mergeTrials = 10;

                    Mr2 = zeros(nT/mergeTrials,nB);

                    for i = 1:mergeTrials:nT

                        meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end),u,:)),1);

                        Mr2(j,:) = meanb;

                        j = j+1;

                    end
                else
                    Mr2=Mr(:,u,:);
                    mergeTrials =1;
                end


                [T,B] = size(Mr2);
                j=1;
                for i = 1:trialsPerCath/mergeTrials:T
                    %Build raster
                    M = Mr2(i:min(i+trialsPerCath/mergeTrials-1, end),:).*(1000/bin);
                    [nTrials,nTimes]=size(M);
                    nexttile
                    imagesc((1:nTimes),1:nTrials,squeeze(M));colormap(flipud(gray(64)));
                    xline(preBase/bin, LineWidth=1.5, Color="#77AC30");
                    xline((stimDur+preBase)/bin, LineWidth=1.5, Color="#0072BD");
                    xticks([preBase/bin (round(stimDur/100)*100+preBase)/bin]);
                    xticklabels(xticks*bin)
                    if nSize >1
                        yline([trialsPerCath/mergeTrials/nSize:trialsPerCath/mergeTrials/nSize:trialsPerCath/mergeTrials-1]+0.5,LineWidth=1)
                    end
                    %caxis([0 max(Mr2,[],'all').*(1000/bin)]);
                    %caxis([0 1]);
                    set(gca,'YTickLabel',[]);

                    if i < T - (trialsPerCath/mergeTrials)*max(positionsMatrix(:))-1
                        set(gca,'XTickLabel',[]);

                    end

                    %
                    %                     if posX(u,j) == NeuronVals(u,respVali(u),3) & posY(u,j) == NeuronVals(u,respVali(u),2)
                    %                         rectangle('Position', [round(posX(u,j)/bin)+round(preBase/bin), 0, round(duration/bin), trialsPerCath/mergeTrials+1],...
                    %                             'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
                    %                     else
                    %                         rectangle('Position', [round(posX(u,j)/bin)+round(preBase/bin), 0, round(duration/bin), trialsPerCath/mergeTrials+1],...
                    %                             'EdgeColor', 'b', 'LineWidth', 1,'LineStyle','-.');
                    %                     end

                    offsetR = 300;
%                     if posX(u,j) == NeuronVals(u,respVali(u),3) & posY(u,j) == NeuronVals(u,respVali(u),2)
%                         rectangle('Position', [round(offsetR/bin)+round(preBase/bin), 0, round(duration/bin), trialsPerCath/mergeTrials+1],...
%                             'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
%                     else
% %                         rectangle('Position', [round(offsetR/bin)+round(preBase/bin), 0, round(duration/bin), trialsPerCath/mergeTrials+1],...
% %                             'EdgeColor', 'b', 'LineWidth', 1,'LineStyle','-.');
%                     end


                    j = j+1;
                end
                fig = gcf;
                set(fig, 'Color', 'w');
%                 c = colorbar;
% 
%                 ylabel(c, 'spikes/sec','FontSize',10);
                %c.Ticks = [0 round((max(Mr2,[],'all')/2)) round(max(Mr2,[],'all'))];
                %c.TickLabels = [0 round((max(Mr2,[],'all')/2)*1000/bin) round(max(Mr2,[],'all')*1000/bin)];

                %%Testing max window selec
% 
%                 figure;imagesc(squeeze(Mr2))
%                 hold on
%                 for j = 1:81
%                 
%                     if posX(u,j) == NeuronVals(u,respVali(u),3) & posY(u,j) == NeuronVals(u,respVali(u),2)
%                         rectangle('Position', [round(posX(u,j)/bin)+round(preBase/bin), j*2, round(duration/bin), 2],...
%                             'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
%                     else
%                         rectangle('Position', [round(posX(u,j)/bin)+round(preBase/bin), j*2, round(duration/bin), 2],...
%                             'EdgeColor', 'b', 'LineWidth', 1,'LineStyle','-.');
%                     end
% %%%Take mean every category
% 
%                 end
                % Set the color of the figure and axes to black
                colorbar;
                title(t,sprintf('Rect-GRid-raster-U%d',u))
                ylabel(t,sprintf('%d trials',nTrials*mergeTrials))
                xlabel(t,'Time (ms)')
                fig.Position =  [147   270   662   446];%[147    58   994   658];
                %prettify_plot
                print(fig,sprintf('%s-rect-GRid-raster-U%d.png',NP.recordingName,u),'-dpng')
                close

            end

        end
    end
%
    if heatMap ==1
    % Build heat map

    cd(NP.recordingDir)
    if ~exist(path+"\Figs",'dir')
        mkdir Figs
    end
    cd(NP.recordingDir + "\Figs")

    trialDiv  = length(seqMatrix)/length(unique(seqMatrix))/nSize;

    offsetR=50;
    duration = 300;%ms for offset response

    bin =1;
    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted+offsetR)/bin),round(duration)/bin);
    [Mro] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted+stimDur)/bin),round(duration)/bin);
    [Mb1] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-duration)/bin),round(duration)/bin);

    [nT,nN,NB] = size(Mr);
    [nTo,nNo,NBo] = size(Mro);
   

    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin));

    Nb2=  mean(Mb,3);

    Nbase = mean(Mb,[1 3]);


    %%%%%%%%%%%%%%%%%%%% Shuffle raster before convolution in order
    %%%%%%%%%%%%%%%%%%%% to calculate tuning index
    %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

%     nShuffle =50;
% 
%     shuffledDataOn = zeros(size(MrMean,1), size(MrMean,2), size(MrC,3), nShuffle);
%     shuffledDataOff = zeros(size(Mro,1), size(Mro,2), size(Mro,3), nShuffle);
% 
%     for i =1:nShuffle
% 
%         % Shuffle along the first dimension
%         idx1 = randperm(size(Mro,1));
% 
%         % Shuffle along the third dimension
%         idx3 = randperm(size(Mro,3));
% 
%         shuffledDataOff(:,:,:,i) = Mro(idx1, :, idx3);
%         shuffledDataOn(:,:,:,i) = Mr(idx1, :, idx3);
% 
%     end


    if noEyeMoves
        % EyePositionAnalysis
        % Create spike Sums with NaNs when the eye is not present.


        file = dir (NP.recordingDir);
        filenames = {file.name};
        files= filenames(contains(filenames,"timeSnipsNoMov-31"));
        cd(NP.recordingDir)
        %Run eyePosition Analysis to find no movement timeSnips
        timeSnips = load(files{1}).timeSnips;
        timeSnipsMode = timeSnips(:,timeSnips(3,:) == mode(timeSnips(3,:)));

        selecInd = [];
        for i = 1:size(timeSnipsMode,2)

            %Find stimOns and offs that are between each timeSnip
            selecInd = [selecInd find(directimesSorted>=timeSnipsMode(1,i) & directimesSorted<(timeSnipsMode(2,i)-stimDur))];
        end

        %MrC = nan(round(nT/trialDiv),nN, NB+NBo);

        MrC = nan(round(nT/trialDiv),nN, NB);


        %%Create summary of identical trials

        MrMean = nan(round(nT/trialDiv),nN);

        for u = 1:length(goodU)
            j=1;

            for i = 1:trialDiv:nT

                indexVal = selecInd(selecInd>=i & selecInd<=i+trialDiv-1);

                if ~isempty(indexVal)


                    meanRon =  reshape(mean(Mr(indexVal,u,:),1),[1,size(Mr,3)]); 

                    meanRoff =  reshape(mean(Mro(indexVal,u,:),1),[1,size(Mro,3)]);

                    %meanBase =  reshape(mean(Mb1(indexVal,u,:),1),[1,size(Mb1,3)]);

                    %MrC(j,u,:) = [meanRon-meanBase meanRoff-meanBase]; %Combine on and off response and substract to each the mean baseline

                    MrC(j,u,:) = [meanRon];
                    MrMean(j,u) = mean(MrC(j,u,:),3);%-Nbase;

                else
                    2+2
                end

            
                j = j+1;

            end
        end

        2+2


    else


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
        
        MrMean = mean(MrC,3);%-Nbase;

    end

  

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

        if string(rect.VSMetaData.allPropVal{find(strcmp(rect.VSMetaData.allPropName,'shape'))}) == "circle"  %%%Check
            
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

            %figure;imagesc(xyScreen')

            VideoScreen(:,:,j) = xyScreen';
            
        elseif nSize ==1 %%Add the possibility of several sizes but square.
            try %New format
                X2 = round(rect.VSMetaData.allPropVal{prop,1}.X2{1}(C(i,2)));
                X1 = round(rect.VSMetaData.allPropVal{prop,1}.X1{1}(C(i,2)));
                X3 = round(rect.VSMetaData.allPropVal{prop,1}.X3{1}(C(i,2)));
                X4 = round(rect.VSMetaData.allPropVal{prop,1}.X4{1}(C(i,2)));
                X = [X1,X2,X3,X4]./reduF;

                Y4 = round(rect.VSMetaData.allPropVal{prop,1}.Y4{1}(C(i,2)));
                Y1 = round(rect.VSMetaData.allPropVal{prop,1}.Y1{1}(C(i,2)));
                Y3 = round(rect.VSMetaData.allPropVal{prop,1}.Y3{1}(C(i,2)));
                Y2 = round(rect.VSMetaData.allPropVal{prop,1}.Y2{1}(C(i,2)));
                Y = [Y1,Y2,Y3,Y4]./reduF;
            catch %old format
                X2 = round(rect.VSMetaData.allPropVal{prop,1}.X2(C(i,2)));
                X1 = round(rect.VSMetaData.allPropVal{prop,1}.X1(C(i,2)));
                X3 = round(rect.VSMetaData.allPropVal{prop,1}.X3(C(i,2)));
                X4 = round(rect.VSMetaData.allPropVal{prop,1}.X4(C(i,2)));
                X = [X1,X2,X3,X4]./reduF;

                Y4 = round(rect.VSMetaData.allPropVal{prop,1}.Y4(C(i,2)));
                Y1 = round(rect.VSMetaData.allPropVal{prop,1}.Y1(C(i,2)));
                Y3 = round(rect.VSMetaData.allPropVal{prop,1}.Y3(C(i,2)));
                Y2 = round(rect.VSMetaData.allPropVal{prop,1}.Y2(C(i,2)));
                Y = [Y1,Y2,Y3,Y4]./reduF;

            end

            mask = poly2mask(X, Y, screenRed, screenRed)';
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

    NanPos = isnan(MrMean);

    MrMean(NanPos) = 0;

    Res = reshape(MrMean,[1,1,size(MrMean,1),nN]).*1000;

    RFu = squeeze(sum(VD.*Res,3));
 

    %%%Calculate Z-score of responses 
    N_bootstrap = 1000;
    cd(NP.recordingDir)

    %boot_means = load(sprintf('RectGrid-Base-Boot-1000-%s',NP.recordingName)).boot_means;
    normMatrixMean = repmat(pxyScreen,[1,1,nN]).*reshape(Nbase,[1,1,nN]);
    normMatrixSTD = repmat(pxyScreen,[1,1,nN]).*((reshape(std(Nb2),[1,1,nN]))+eps);%&+1/(nT));

    %normRFu = (RFu-repmat(pxyScreen,[1,1,nN]).*mean)./normMatrix;

    normRFu = (RFu-normMatrixMean)./normMatrixSTD;
    %normRFu = (RFu)./repmat(pxyScreen,[1,1,nN]);

    %figure;imagesc(RFu(:,:,60))

    save(sprintf('RFuStatic-%s',NP.recordingName),"RFu")
%&

    if plotHeatMap
         eNeuron =[72];%selecN{1}(2,selecN{1}(1,:)==ex);
         cd(NP.recordingDir)
         %Parameters
         eye_to_monitor_distance = 21.5; % Distance from eye to monitor in cm
         pixel_size = 33/(1080/reduF); % Size of one pixel in cm (e.g., 25 micrometers)
         monitor_resolution = [screenRed, screenRed]; % Width and height in pixels
         [theta_x theta_y] = pixels2eyeDegrees(eye_to_monitor_distance,pixel_size,monitor_resolution);

         if noEyeMoves
            ResNanPos = reshape(NanPos,[offsetN, offsetN, size(NanPos,2)]);
         else
             zeros = size(offsetN,offsetN,nN);
         end

         
         for u = eNeuron
             ResNanPosU = ResNanPos(:,:,u)';
             M = normRFu(:,:,u);
             M(ResNanPosU) = 0;
             fig = figure;imagesc(M)
             cd(NP.recordingDir)
             colorbarMBlims = load(sprintf('%s-Unit-%d-MovBall-RFlims-Dirs.mat',NP.recordingName,u)).colorbarlims;
             c = colorbar;title(c,'Z-score');%caxis([0 6.2])
             set(gcf,'Color','w')
             colormap('jet')
             xt = xticks;
             xticklabels(round(theta_x(1,xt)))
             yt = yticks;
             yticklabels(round(theta_y(yt,1)))
             xlabel('X degrees')
             ylabel('Y degrees')
             fig.Position = [ 527   511   455   360];
             cd(NP.recordingDir+"\Figs")
             print(gcf,sprintf('%s-NEM-Unit-%d-rectGrid-RF.png',NP.recordingName,u),'-dpng')


         end

    end
%
    if calculateEntropy ==1

        offsetN = cell2mat(rect.VSMetaData.allPropVal(find(strcmp(rect.VSMetaData.allPropName,'rectGridSize'))));
        offsetN = 9; %Make it the same for all!!

         entropies = zeros(1,nN);

         if noEyeMoves
             ResNanPos = reshape(NanPos,[offsetN, offsetN, size(NanPos,2)]);
         else
             ResNanPos = zeros(offsetN,offsetN,nN);
         end

            for u = 1:nN
                
                M = squeeze(normRFu(:,:,u));
                

                ResNanPosU = ResNanPos(:,:,u)';

                %
                %figure;imagesc(P_scaled);colorbar

                %%Prepare matrix for entropy calculation


                % Define edges for rows and columns
                rowEdges = round(linspace(1, size(M, 1) + 1, offsetN + 1));
                colEdges = round(linspace(1, size(M, 2) + 1, offsetN + 1));

                % Initialize the resulting matrix
                reducedMatrix = zeros(offsetN, offsetN);

                % Assign central element of each block
                for i = 1:offsetN
                    for j = 1:offsetN
                        % Extract block using edges
                        block = M(rowEdges(i):rowEdges(i+1)-1, colEdges(j):colEdges(j+1)-1);

%                        % Find the central element of the block
%                         centralRow = ceil(size(block, 1) / 2);
%                         centralCol = ceil(size(block, 2) / 2);
%                         centralValue = block(centralRow, centralCol);

                        % Assign the central value to the reduced matrix
                        reducedMatrix(i, j) = mean(block(:), 'omitnan');
                    end
                end


                if noEyeMoves
                    reducedMatrix(ResNanPosU) = 0;
                end
                
                
              
                % Normalize to create a probability distribution
               % reducedMatrix = reducedMatrix / sum(reducedMatrix(:),'omitnan');

                % Convert to an image-like format and calculate entropy
                % (scale to [0, 1] for compatibility with `entropy`)
                P_scaled = mat2gray(reducedMatrix);
                %figure;imagesc(reducedMatrix);colorbar;
                entropies(u) = entropy(P_scaled);
            end

            cd(NP.recordingDir)
            %sign = 0.005;
            if noEyeMoves
                save(sprintf('NEM-Entropies-RG-RF-respU-%s.mat',NP.recordingName),'entropies')
            else
                save(sprintf('Entropies-RG-RF-respU-%s.mat',NP.recordingName),'entropies')
            end



    end




    end
end












