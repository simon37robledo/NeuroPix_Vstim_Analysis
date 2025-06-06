cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%Optional
bombcelled = 0;
sigma = 0.05; %Criteria for responsive units. Used on the tunin and receptive field analysis. 

summPlot = 0;
plotexamplesMBR =1;
newTIC = 0;

ResponseProfile=1; redoResp=0;
Shuffling =0;
Shuffling_baseline=0;%Everything that involves the TIC matrix needs to change (choose trials) 
repeatShuff =0;
trialThres =0.6;

ReceptiveFieldFixedDelay = 0;

tuning =0;
ZscoreMethod = 0;


depthPlot =0;

ReceptiveFieldConvolutions =1;
repeatConv =1;
useZscore = 0;
noEyeMoves = 1;
ModeQuadrant = 1;
XYdivision =0;
plotRF =0;

spatialTuning=0;
calculateEntropy =0;

includeOnespeed = 1;
x=1;
examplesSDG = [1 2 3 4 5 6 7 29 30 31 32 40 41 42 43];
N_bootstrap = 1000;


%examplesSDG =[1 2 3 4 5 6 7 29 30 31 32 40 41 42 43];

pv27 = [8 9 10 11 12 13 14];

newDiode =0;
GoodRecordingsRF = [8:20,40:43,49:54];
Awake = [44:48];

duration = 300; %in ms. Duration of window in which the response is calculated. This window
% moves along the raster and the maximal mean is selected as response, per trial, offset, etx
%GoodRecordingsPV =[1:21,40:43,49:54];
%Plot specific neurons
 
%%%In shuffling make sure that response cat is selected equally between SDG
%%%and MB
%%
%r=1;%check rect
% Iterate through experiments (insertions and animals) in excel file
for ex =  [51]%GoodRecordingsPV%GoodRecordingsPV%SDGrecordingsA%GoodRecordings%GoodRecordingsPV%GoodRecordingsPV%selecN{1}(1,:) %1:size(data,1)
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

    %Create Figs and matData folders if they don't exist


    if ~exist(path+"\Figs",'dir')
        mkdir Figs
    end

    if ~exist(path+"\matData",'dir')
        mkdir matData
    end


    %2. Extract moving ball statistics
    patternIndex = strfind(string(NP.recordingDir), "\catgt");

    endIndex = patternIndex(1)-1;
    stimDir = string(NP.recordingDir);
    stimDir = extractBetween(stimDir,1,endIndex);

    file = dir (stimDir);
    filenames = {file.name};

    file = dir (stimDir);
    filenames = {file.name};
    ballFiles = filenames(contains(filenames,"linearlyMovingBar"));

    
        if isempty(ballFiles)
            %disp()
            w= sprintf('No moving bar files where found in %s. Skipping into next experiment.',NP.recordingName);
            warning(w)
            %continue
        end

    directions = [];
    lengthBar = [];
    widthBar = [];
    speeds = [];

    j =1;

    if size(ballFiles) ~= [0 0]

        for i = ballFiles %Extract stim properties
            ball= load(stimDir+"\"+string(i));

            directions = [directions cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'directions'))))];

            widthBar = [widthBar cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'barWidth'))))];

            lengthBar = [lengthBar cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'barLength'))))];

            speeds = [speeds cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'speeds'))))];

            interStimStats = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    
    uDir = unique(directions);
    direcN = length(uDir);
    uSpeed = unique(speeds);
    speedN = length(uSpeed);
    %offsetN = length(unique(offsets));
    nT = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials'))));
    trialDivision = nT/(length(uDir)*length(uSpeed)); %Number of trials per unique conditions
    categ = nT/trialDivision;


    %3. Load Triggers (diode)
    Ordered_stims= strsplit(data.VS_ordered{ex},',');
    %containsMB = cellfun(@(x) contains(x,'MB'),Ordered_stims);
    containsMB = strcmp(Ordered_stims, 'MBR');
    ttlInd = find(containsMB);



    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,newDiode,1,"MBR",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,0,1,"MBR",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A


     %When dealing with different speeds, save different stim durs, and create M for each speed
    A = [stimOn directions' speeds'];
    [C indexS] = sortrows(A,[2 3]);

    B = [stimOff directions' speeds'];
    [Coff indexSo] = sortrows(B,[2 3]);

    stimInter= mean(stimOn(2:end)-stimOff(1:end-1));

    if includeOnespeed

        A =  A(A(:,3) == max(A(:,3)),:);

        B =  B(B(:,3) == max(B(:,3)),:);

        [C indexS] = sortrows(A,[2 3]);

        [Coff indexSo] = sortrows(B,[2 3]);

        stimOn = A(:,1);
        stimOff = B(:,1);

        %digon  =onDigital(A(:,5) == max(A(:,5)));

        stimInter = interStimStats;

    end

    %4. Sort directions:
    directimesSorted = C(:,1)';
    directimesSortedOff = Coff(:,1)';

    % When dealing with different speeds, save different stim durs, and create M for each speed
    stimDurs = -stimOn+stimOff;

    stimDur = mean(stimDurs(stimDurs>mean(stimDurs)-std(stimDurs)));

    %5. Load data to select good units
    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');
    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");
    GoodU_orDepth = cluster_info.depth(cluster_info.group=="good");

    %If new tic matrix needs to be used, move oldTIc files and run convertPhySorting2tIc
    cd(NP.recordingDir)
    if isfile("sorting_tIc.mat") && newTIC

        if ~exist('oldTIC', 'dir')
            mkdir oldTIC
        end

        movefile sorting_tIc.mat oldTIC
        movefile sorting_tIc_All.mat oldTIC
    end

    if bombcelled
        p = NP.convertPhySorting2tIc(NP.recordingDir,0,1);
    else
        p = NP.convertPhySorting2tIc(NP.recordingDir);
    end

   
    
    label = string(p.label');
    goodU = p.ic(:,label == 'good');
    phy_IDg = p.phy_ID(label == 'good');
    
    if isempty(goodU)
        %disp()
        w= sprintf('No somatic neurons in %s. Skipping into next experiment.',NP.recordingName);
        warning(w)
       % continue
    end
    goodSPKS = p.nSpks(label == 'good');
    goodAmp = p.neuronAmp(label == 'good');
    goodUdepth = NP.chLayoutPositions(2,goodU(1,:));

    %7. Load raster matrices
    bin = 1;
    preBase = round(stimInter-200);%round(3*interStimStats/4);

     if noEyeMoves
% 
%         %%%%%% Construct stimType matrix for eye movement plot.
%         stimType = zeros(length(A),6); %3 plus number of example neurons
%         stimType(:,1) = A(:,1);
%         stimType(:,2) = A(:,1)+stimDur;
%         stimType(:,3) = A(:,2);
%         stimType(:,4) = A(:,3);
% 
%         %Get response strenght of specific neurons and save it in stimType
%         [MrNoSort] = BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn'-preBase)/bin),round((stimDur+preBase*2)/bin)); %response matrix
%         ResponseStrengthU34 = mean(MrNoSort(:,34,round(preBase/bin):round((preBase+stimDur)/bin)),3); %For PV35_3
%         ResponseStrengthU8 =  mean(MrNoSort(:,8,round(preBase/bin):round((preBase+stimDur)/bin)),3); %For PV35_3
%         stimType(:,end-1) = ResponseStrengthU34;
%         stimType(:,end) = ResponseStrengthU8;
% 
           EyePositionAnalysis(NP,data.Eye_video_dir{ex},31)
% 
     end
    

    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-stimInter/2)/bin),round((stimDur+stimInter)/bin)); %response matrix
    [nT,nN,nB] = size(Mr);

    %figure;imagesc(squeeze(Mr(:,15,:)));colormap(flipud(gray(64)));xline(preBase/bin);xline((preBase+stimDur)/bin);yline(trialDivision:trialDivision:nT-1);

    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((preBase)/bin));% baseline matrix (whole baseline)
    baseline = squeeze(mean(squeeze(mean(Mb)),2));
    MrMean = squeeze(mean(Mr));

    v_replicated = repmat(baseline, 1, size(MrMean, 2));
    Norm = MrMean ./ v_replicated; %normalization by division
    Norm(Norm>2) = 2;

    %8. Plot summary response
%
    for plotOp = summPlot
        if summPlot

            for s = 1:speedN
                S = Coff(:,3)'== uSpeed(s);
                stimDurS = mean(-directimesSorted(S)+directimesSortedOff(S));

                [MrS] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted(S)-stimInter/2)/bin),round((stimDurS+stimInter)/bin)); %response matrix
                [nT,nN,nB] = size(MrS);

                %figure;imagesc(squeeze(Mr(:,15,:)));colormap(flipud(gray(64)));xline(preBase/bin);xline((preBase+stimDur)/bin);yline(trialDivision:trialDivision:nT-1);

                [MbS] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted((S))-preBase)/bin),round((preBase)/bin));% baseline matrix (whole baseline)
                
                baseline = squeeze(mean(squeeze(mean(MbS)),2));
                MrMean = squeeze(mean(MrS));

                v_replicated = repmat(baseline, 1, size(MrMean, 2));
                Norm = MrMean ./ v_replicated; %normalization by division
                Norm(Norm>5) = 5;

                fig = figure;
                imagesc(Norm);
                xline(preBase/bin,'k', LineWidth=1.5)
                xline(stimDurS/bin+preBase/bin,'k',LineWidth=1.5)
                hcb = colorbar();
                title(hcb,'SpkR/Baseline');
%                 xticks([0.5 (preBase/bin):50:nB])
%                 xticklabels([-preBase 0:10*bin:nB*bin])
                xt = xticks;
                xticklabels(round(xticks/1000))
                ylabel('Neurons');xlabel('Time (ms)');
                % Define key colors: blue, white, yellow
                keyColors = [0 0 0.5; 1 1 1;[0 0.5 0]]; % RGB for blue, white, yellow
                % Number of colors to interpolate between each key color
                nInterpolations = 32;
                % Interpolate colors between blue and white
                blueToWhite = interp1([1,2], keyColors(1:2,:), linspace(1,2,nInterpolations));
                % Interpolate colors between white and yellow
                whiteToYellow = interp1([1,2], keyColors(2:3,:), linspace(1,2,nInterpolations));
                % Combine the two gradients, omitting one instance of white to maintain 64 colors
                customColormap = [blueToWhite; whiteToYellow(2:end,:)];
                % Apply the custom colormap
                colormap(customColormap);
                title(sprintf('Speed = %d',uSpeed(s)))
                cd(NP.recordingDir)
                cd(NP.recordingDir + "\Figs")
                print(fig, sprintf('%s-MovBall-summary.png',NP.recordingName),'-dpng');
            end

        end
    end

    stims = stimOn';

    %%%%%%%%%%%%%% Select baseline
    %NP.recordingDuration_ms/1000/60

    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

    Mb = mean(Mb,3); %mean across time bins

% 
     spkRateBM = mean(Mb); %Calculate baseline spike rate of each neuron.

    %2. Calculate denominator for Z-score calculation
    %spkRateBM = mean(Mb); %total mean.
    epsilon = 0.01;
    %denom = mad(MbC,0)+epsilon; %Calculate baseline variation of each neuron.

    denom = mad(Mb,0)+eps;
    %%%%%%%%%%%%%%%%%%%% Select responsive units and calculate their tunning profile
    %3. Convolute matrix to spread response
    % Convolute in the 3rd dimension (trials)
    for ResponseProfileF =1
if ResponseProfile ==1
    %4. Create window to scan rasters and get the maximum response
    
    Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((stimDur+duration)/bin)); %response matrix
    [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');

    [nT,nN,nB] = size(MrC);

   
    window_size = [1, round(duration/bin)];

    %5. Initialize the maximum mean value and its position

    max_position = zeros(nN,2);
    max_mean_value = zeros(1,nN);
    max_mean_valueB = zeros(1,nN);

    NeuronVals = zeros(nN,nT/trialDivision,6); %Each neuron, which has a matrix where the first column is maxVal of bins, 2nd and 3rd position of window in matrix...
    % 4th Z-score.
    % responsive compared to the baseline.

    %[1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 17, 20, 22, 23, 24, 28, 32, 34, 36]

    if includeOnespeed
        responseFile = sprintf('OneSpeed-%d-NeuronRespCatBar-%s.mat',max(speeds),NP.recordingName);
    else
        responseFile = sprintf('NeuronRespCatBar-%s.mat',NP.recordingName);
    end

    if ~isfile(sprintf('NeuronRespCatBar-%s.mat',NP.recordingName)) || redoResp
    
    %Baseline = size window

    [Mbd] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

    %Merge trials:

    mergeTrials = trialDivision;

    Bd = reshape(Mbd, [mergeTrials, size(Mbd, 1)/mergeTrials, size(Mbd, 2), size(Mbd,3)]);

    Mbd2 = squeeze(mean(Bd, 1));


    B = reshape(MrC, [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);

    % Take the mean across the first dimension (rows)
    Mr2 = squeeze(mean(B, 1));

    [nT,nN,nB] = size(Mr2);
    cd(NP.recordingDir)

    %Real data:
    
    for u =1:nN
        % Slide the window over the matrix
        %unit matrix
        max_mean_value(u) = -Inf; %General max? not needed
        max_mean_valueB(u)=-Inf;
        NeuronRespProfile = zeros(nT,6); %4 columns plus: , dir, , speed, .

        k =1;
        max_position_Trial = zeros(nT,2); %Create 2 matrices, for mean inside max window, and for position of window, for each trial category
        max_mean_value_Trial = zeros(1,nT);
        max_mean_value_TrialB = zeros(1,nT);
        for i = 1:nT %Iterate across trials %%%%%%%%%% 
            uM = squeeze(Mr2(i,u,:))';%Create matrix per unique trial conditions

            uMB = squeeze(Mbd2(i,u,:))';%Create matrix per unique trial conditions
            %uMb = Mbd2(i,u);

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
            %Assign visual stats to last columns of NeuronRespProfile. Select
            %according  to trial (d) the appropiate parameter > directions'offsets' sizes' speeds' freq'
            NeuronRespProfile(k,5) = C(i*mergeTrials,2);
            NeuronRespProfile(k,6) = C(i*mergeTrials,3);

            k = k+1;

        end

        %figure;imagesc(uM);xline(max_position_Trial(i,2));xline(max_position_Trial(i,2)+window_size(2))

        NeuronVals(u,:,:) = NeuronRespProfile;
    end
            if includeOnespeed && length(unique(speeds)) > 1
                save(sprintf('OneSpeed-%d-NeuronRespCatBar-%s',max(speeds),NP.recordingName),"NeuronVals")
            else
                save(sprintf('NeuronRespCatBar-%s',NP.recordingName),"NeuronVals")
            end
    else

        if includeOnespeed && length(unique(speeds)) > 1
            NeuronVals = load(sprintf('OneSpeed-%d-NeuronRespCatBar-%s',max(speeds),NP.recordingName)).NeuronVals;
        else
            NeuronVals = load(sprintf('NeuronRespCatBar-%s',NP.recordingName)).NeuronVals;
        end
    end

    %Create tunning curve, based on direction tunning curve 
    tuningCurve = zeros(nN,direcN);


    udir = unique(directions);

    for u = 1:nN
        NeuronD = squeeze(NeuronVals(u,:,[1 5])); %1, mean response, 5, direction.
        for d = 1:direcN
            tuningCurve(u,d) = max(NeuronD(NeuronD(:,2)==udir(d),1))'; %Selecting top direction. 
        end
    end

    if ~isfile(sprintf('NeuronRespCat-%s.mat',NP.recordingName)) || redoResp

        save(sprintf('tuningC-MBR-%s',NP.recordingName),"tuningCurve")

    else

        tuningCurve = load(sprintf('tuningC-MBR-%s',NP.recordingName)).tuningCurve;

    end

    
for Shuffle =1
    if Shuffling ==1

        % %Keep real max window for shuffling.
        % [mp mi] = max(NeuronVals(:,:,1),[],2);
        %
        % psT = zeros(size(mi, 1), 1);
        % psB = zeros(size(mi, 1), 1);
        %
        % % Extract the values from A(:,:,2:3) according to the max indices
        % for i = 1:size(mi, 1)
        %     psT(i) = NeuronVals(i, mi(i), 2);
        %     psB(i) = NeuronVals(i, mi(i), 3);
        % end
        rands = 1000;
        [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-stimInter/2)/bin),...
            round((stimDur+stimInter)/bin)); %response matrix
        sMr = single(Mr);
        %sMb = single(Mbd);
%         [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round(preBase/bin));
%         sMb = single(Mb);
        if ~isfile(sprintf('randValues-It-%d.mat',rands))||repeatShuff==1
            %Select significantly responsive units with shuffling plus Z-score:
           
            %            numSegments = 4;
            %             parfor par = 1:numSegments
            RandValU = zeros(rands, nN,2,'single');
            tic

            for s = 1:rands
                    tic

                    %%%Shuffle trials
                    %check if neuron is responsive to specific stim (shuffle across trials)
                    %M_shufMTrC = ConvBurstMatrix(M_shuffTr,fspecial('gaussian',[1 5],3),'same');
                    B = reshape(sMr(randperm(size(Mr,1)),:,:), [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);
                    M_shuffTr =  squeeze(mean(B, 1));

                    %%%Shuffle times
                    %             ti2 = zeros(1,length(p.t));
                    %             for u = 1:nN
                    %                 ti=p.t(p.ic(3,u):p.ic(4,u));
                    %                 diffi = diff(ti);
                    %                 ti = ti(1) + cumsum(randperm(round(diff(ti))));
                    %                 ti2(p.ic(3,u):p.ic(4,u)) = [ti(1) ti(1)+cumsum((diffi(randperm(numel(diffi)))))];
                    %             end
                    %
                    %             M_shuffTi=BuildBurstMatrix(goodU,round(ti2/bin),round((stims-preBase)/bin),round((stimDur+preBase*2)/bin)); %check if neuron is visually response (shuffle across time)
                    %             M_shufMTiC = ConvBurstMatrix(M_shuffTi,fspecial('gaussian',[1 5],3),'same');

                    %B = reshape(sMr(randperm(size(Mr,1)),:,randperm(size(Mr,3))), [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);
                    B = reshape(sMr(randperm(size(Mr,1)),:,randperm(size(Mr,3))), [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);
                    M_shuffTi =  squeeze(mean(B, 1));

                    %                 %MB_shufTi = BuildBurstMatrix(goodU,round(ti2/bin),round((stims-preBase)/bin),round(preBase/bin));
                    %                 %MB_shufTi = mean(MB_shufTi,3);
                    %                 spkRateBMS = mean(MB_shufTi); %total mean.
                    %                 denom = mad(MbC,0)+epsilon; %Calculate baseline variation of each neuron.

                    %denomSti = std(MB_shufTi)+eps;
                    
                    
                    NeuronRespProfileShTr = zeros(nN,'single');
                    NeuronRespProfileShTi = zeros(nN,'single');
                    NeuronRespProfileShDiff = zeros(nN,'single');


                    for u =1:nN
                        %                     Slide the window over the matrix (this selects the random
                        %                     max)
                        %                     unit matrix

                        k =1;
                        for i = 1:nT/trialDivision %Iterate across trials
                            uMTr = squeeze(M_shuffTr(i,u,:))';%Create matrix per unique trial conditions
                            uMTi = squeeze(M_shuffTi(i,u,:))';%Create matrix per unique trial conditions
                            %uMB = Mbd2(i,u);

                            %Create 2 matrices, for mean inside max window, and for position of window, for each trial
                            max_mean_value_Trialtr = zeros(1,nT);
                            max_mean_value_Trialtr(i) = -Inf;

                            %Create 2 matrices, for mean inside max window, and for position of window, for each trial
                            max_mean_value_Trialti = zeros(1,nT);
                            max_mean_value_Trialti(i) = -Inf;

                            for j = 1:2:size(uM, 2) - window_size(2) + 1 %Iterate across bins
                                % Extract the sub-matrix
                                sub_matrixTr = uMTr(j:min(j+window_size(2)-1,end)); %Create matrix of size window per bin
                                sub_matrixTi = uMTi(j:min(j+window_size(2)-1,end)); %Create matrix of size window per bin

                                % Compute the mean value
                                mean_valueTr = mean(sub_matrixTr(:)); %Compute mean of each window
                                mean_valueTi = mean(sub_matrixTi(:)); %Compute mean of each window

                                % Update the maximum mean value and its position (a
                                % window is selected across each trial)
                                if mean_valueTr >  max_mean_value_Trialtr(i)
                                    max_mean_value_Trialtr(i) = mean_valueTr;
                                end

                                % Update the maximum mean value and its position (a
                                % window is selected across each trial)
                                if mean_valueTi >  max_mean_value_Trialti(i)
                                    max_mean_value_Trialti(i) = mean_valueTi;
                                end

                            end

                            %                     NeuronRespProfileShTr(i,u) = (max_mean_value_Trialtr(i) - spkRateBM(u))/denom(u);%Mb_shufTrM(u))/denomTr(u); %Zscore
                            %                     NeuronRespProfileShTi(i,u) = (max_mean_value_Trialti(i) - spkRateBMS(u))/denomSti(u);%Mb_shufTiM(u))/denomTi(u); %Zscore
                            NeuronRespProfileShTr(i,u) = max_mean_value_Trialtr(i);
                            NeuronRespProfileShTi(i,u) = max_mean_value_Trialti(i);
                            NeuronRespProfileShDiff(i,u) = max_mean_value_Trialti(i)-spkRateBM(u);
                            k=k+1;

                            %                     NeuronRespProfileShTr = squeeze(mean(M_shuffTr(psT(u),u,psB(u):psB(u)+round(duration/bin))));
                            %                     NeuronRespProfileShTi = squeeze(mean(M_shuffTi(psT(u),u,psB(u):psB(u)+round(duration/bin))));
                        end

                    end

                    %             tDivList = 1:trialDivision:nT-1;
                    %
                    %             meanTshuffledTr = zeros(u,numel(tDivList));
                    %             meanTshuffledTi = zeros(u,numel(tDivList));
                    %
                    %             %
                    %             for td = 1:numel(tDivList) %take means across trial divisions to look at random times
                    %                 meanTshuffledTi(:,td) = mean(NeuronRespProfileShTi(tDivList(td):tDivList(td)+trialDivision-1,:));
                    %                 meanTshuffledTr(:,td) = mean(NeuronRespProfileShTr(tDivList(td):tDivList(td)+trialDivision-1,:));
                    %             end

                    %tri = randi(numel(tDivList)); %Take group of random trials to look at trials

                    %RandZscore matrix: 3 dimension, 1st (number of iterations) 2nd
                    %(number of neurons), 3rd (1. trials shuffled, 2. Times
                    %shuffled).

                    RandValU(s,:,1) =  max(NeuronRespProfileShTr);
                    RandValU(s,:,2) =  max(NeuronRespProfileShTi);
                    RandValU(s,:,3) = max(NeuronRespProfileShDiff);

                    % RandZscoreU(s,:,2) =  mean(NeuronRespProfileShTi(tDivList(tri):tDivList(tri)+trialDivision-1,:));%max(meanTshuffledTi,[],2);
                    toc
            end


            %             end
            toc
            %             RandValU = cat(1, RandValUp{:});

            save(sprintf('randValues-It-%d',rands),"RandValU");
        else %if eand Z-score file exists then load
            RandValU = load(sprintf('randValues-It-%d.mat',rands)).RandValU;
        end
        toc
        Zthreshold = 3;
        maxZScores = squeeze(max(tuningCurve,[],2));
        %%Check if at least for 1 direction neuron is visually responsive:

        %maxZS_vr = max(meanZScores,[],2);

        pvalTr = (sum(squeeze(RandValU(:,:,1)) >= maxZScores') +1)/(rands+1);
        pvalTi = (sum(squeeze(RandValU(:,:,2)) >= maxZScores') +1)/(rands+1);
        pvalDiff = (sum(squeeze(RandValU(:,:,3)) >= maxZScores') +1)/(rands+1);
        %figure;histogram(squeeze(RandZscoreU(:,15,1)))

        save(sprintf('pvalTrials-%s',NP.recordingName),'pvalTr')
        save(sprintf('pvalTime-%s',NP.recordingName),'pvalTi')
        save(sprintf('pvalDiff-%s',NP.recordingName),'pvalDiff')

 
    end

    %%%%%%%%%%%%%%%%%%%%% Botstraping baseline

    if Shuffling_baseline
        
        

        if ~isfile(sprintf('pvalsBaselineBootMBR-%d-%s.mat',N_bootstrap,NP.recordingName))||repeatShuff==1

        baseline = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((preBase)/bin));
        baseline = single(baseline);
        [nT,nN,nB] = size(baseline);
        
        % Bootstrapping settings
      
        boot_means = zeros(N_bootstrap, nN,'single');
        resampled_indicesTr = single(randi(nT, [trialDivision, N_bootstrap]));% To store bootstrapped means

%         temp = NaN(trialDivision,duration);
%         rand_lev1 = randi(num_lev1,num_lev1,1);
%         for j = 1:length(rand_lev1)
%             num_lev2 = find(~isnan(data(rand_lev1(j),:)),1,'last'); %We need to calculate this again here because there is a different number of trials for each neuron
%             rand_lev2 = randi(num_lev2,1,num_trials); %Resample only from trials with data but same number of sample trials for all
%             temp(j,:) = data(rand_lev1(j),rand_lev2);
%         end
%     
        
        resampled_indicesTi = single(randi(nB, [duration, N_bootstrap]));

        kernel = ones(trialDivision, duration) / (trialDivision * duration); % Normalize for mean
        % Start a parallel pool (if not already started)
%         if isempty(gcp('nocreate'))
%             parpool; % Start a pool with the default number of workers
%         end

        tic
        for i = 1:N_bootstrap
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

        %[bootstats] = get_bootstrapped_equalsamples(data,nruns,num_trials,param)

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


            if emptyRows/trialDivision >= trialThres %%%Check which unit is 34 i respU
                pvalsResponse(u) = 1;
            end

        end
        save(sprintf('pvalsBaselineBootMBR-%d-%s',N_bootstrap,NP.recordingName),'pvalsResponse')
        save(sprintf('MovBall-ZscoreBootMBR-%d-%s',N_bootstrap,NP.recordingName),'ZScoreU')
        save(sprintf('MovBall-Base-BootMBR-%d-%s',N_bootstrap,NP.recordingName),'boot_means')

        end

    end

end
end
    
    end %End of Z-score data


%     find(pvalTr<0.05 & pvalTr<0.02)
% 
%     randddd2 = squeeze(RandZscoreU(:,:,1));
% 
%     %Calculate direction: if neuron has 2 big peaks (equivalent angles (hue), or if neuron has 1
%     %big peak (hue), or if neuron has multiple big peaks.
%     binarizedS = meanZScores>Zthreshold;
% 
%     %remove neurons that have zero spiking rate in 97% of the trials
%     rateTrials = mean(Mr,3);
%     tuningCurveS = tuningCurve(sum(rateTrials==0)./size(rateTrials,1)<0.95,:,:);
% 
%     meanZScores = squeeze(mean(tuningCurveS(~isoutlier(tuningCurveS,2)),2));
%     meanZScoresSig =  meanZScores(any(meanZScores > Zthreshold, 2),:);
%     verticalDepthSig = verticalDepth(any(meanZScores > Zthreshold, 2));
% 
%     figure;
%     plot(rad2deg(uDir), meanZScores(1,:), 'LineWidth', 2);
%
%
for tuningIndex = 1
if tuning ==1


    cd(NP.recordingDir)

    pvals= load(sprintf('pvalsBaselineBoot-1000-%s',NP.recordingName)).pvalsResponse;
    pvalsbar= load(sprintf('pvalsBaselineBootMBR-1000-%s',NP.recordingName)).pvalsResponse;
    goodNeurons =  find(pvals <sigma | pvalsbar < sigma);

    trials_for_angle = nT/direcN;



    %%% spike rate tuning curve

    tuningCurve = load(sprintf('tuningC-MBR-%s',NP.recordingName)).tuningCurve;
    tuningCurve = tuningCurve(goodNeurons,:);





    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-stimInter)/bin),round((stimDur+stimInter*2)/bin)); %response matrix
 

    nN = size(tuningCurve,1);
    %%%%% Tuning indexes (OSI, DSI, Prefered angle)
    % Calculate tuning direction and strength for all neurons then plot it across depth and in 3D
    uDirDeg = rad2deg(uDir);
    
    %%% method one OSI %Same result
    OI = sqrt(sum(tuningCurve.*sin(2*uDir),2).^2 + sum(tuningCurve.*cos(2*uDir),2).^2)./sum(tuningCurve,2); % OI calculated as in Dragoi and Swindale.

    %%%% Preferred orientation
    oneSideCurve = zeros(nN,direcN/2);
  

    a = sum(tuningCurve.*cos(2*uDir),2);
    b = sum(tuningCurve.*sin(2*uDir),2);

    %%%Calculate prefered angle
    Theta = zeros(1,nN);
    for u = 1:nN
        if a(u)>0
            Theta(u) = rad2deg(0.5*atan(b(u)/a(u)));
        else
            Theta(u) = rad2deg(deg2rad(180)+0.5*atan(b(u)/a(u)));
        end
    end

    %%% method 2 OSI
    L = sqrt(a.^2+b.^2)./sum(tuningCurve,2); %Tuning Strength

    [preferDir pI] = max(tuningCurve,[],2);

    nonPreferDir = zeros(nN,1);
    for u = 1:nN
        npI = find(uDirDeg == mod(uDirDeg(pI(u)) + 180, 360));
        nonPreferDir(u) = tuningCurve(u,npI);
    end


    DSI = 1- nonPreferDir./preferDir; %%%check rasters of prefer/nprefer =1' 
    DSI(isnan(DSI))=0;
    L(isnan(L))=0;

   % find(DSI(goodNeurons) ==0)

%     if any(DSI == 0) | any(DSI == 0) 
% 2+2
%     end

    Theta = uDirDeg(pI);

    save(sprintf('MBR-Angle-prefer-%s',NP.recordingName),'pI')
    save(sprintf('MBR-Orientation-Tuning-Index-%s',NP.recordingName),'L')
    save(sprintf('MBR-Direction-Selectivity-Index-%s',NP.recordingName),'DSI')
   

    %respU = string(data.ResponseU(ex));




%     % strlength(respU);
%     if strlength(respU)>0
% 
%         strSplit = strsplit(respU, ','); % split the string
%         ResponsiveN = cellfun(@str2num, strSplit); %
%     else
%         ResponsiveN = [];
%     end

    %figure;polarplot(deg2rad(Theta(ResponsiveN)),L(ResponsiveN),'.','MarkerSize',20);set(gcf, 'Color', 'w')


%     find(tuningCurve)



    %figure;imagesc(oneSideCurve)
    %OIn = sqrt(sum(oneSideCurve.*sin(2*uDirRad(1:length(uDir)/2)),2).^2 + sum(oneSideCurve.*cos(2*uDirRad(1:length(uDir)/2)),2).^2)./sum(oneSideCurve,2); % OI calculated as in Dragoi and Swindale.

end

end
%%
% Rasters plots per Neuron
for plotOp = plotexamplesMBR %rstx
    if plotexamplesMBR == 1


        cd(NP.recordingDir)

        posX = squeeze(NeuronVals(:,:,3));
        posY = squeeze(NeuronVals(:,:,2));


        uDir = unique(directions);
        bin = 1;
        bin2 =20;
        trialsPerAngle = trialDivision*speedN;

        %[C sIndex2]= sortrows(A,orderS(k,:));

        %Sort directions:

        %             Co = C;

        %             C = C(C(:,4) == 200,:); %%Select size; PV139_1

        directimesSorted = C(:,1)';

        preBase = round((stimInter));

        [Mr2] = BuildBurstMatrix(goodU,round(p.t/bin2),round((directimesSorted-preBase)/bin2),round((stimDur+preBase*2)/bin2));

        [nT,nN,nB] = size(Mr);


        pvals = load(sprintf('pvalsBaselineBootMBR-%d-%s',N_bootstrap,NP.recordingName),'pvalsResponse').pvalsResponse;
        pvalsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName),'pvalsResponse').pvalsResponse;
        %eNeuron = find(pvals<0.005 |pvalsMB <0.005 );
        eNeuron = find(pvals<0.005);

        sizeN=1;

        for u = eNeuron

            j=1;

            fig =  figure;

            subplot(10,1,[1 6]);
            pos1 = 0.1;
            %subplot('Position', [pos1, 0.1, 0.7, 0.8]);

            imagesc(squeeze(Mr2(:,u,:)).*(1000/bin2));colormap(flipud(gray(64)));
            %Plot stim start:
            xline(preBase/bin2,'k', LineWidth=1.5)
            %Plot stim end:
            xline(stimDur/bin2+preBase/bin2,'k',LineWidth=1.5)
            ylabel(sprintf('%d Trials',nT));
            title(sprintf('U.%d-Unit-phy-%d',u,phy_IDg(u)));
            caxis([0 1])

    
            yline([nT/direcN:nT/direcN:nT-1],'k',LineWidth=3);
 

            yticks([trialDivision:trialDivision:nT])

            yticklabels([repmat([trialDivision:trialDivision:(nT)/4],1,4)])

            hold on
            %Plot rectangle: *max window per unique trial.
            [respVal respVali]= max(NeuronVals(:,:,1),[],2);
            for d = 1:size(NeuronVals,2)

                if posX(u,d) == NeuronVals(u,respVali(u),3) & posY(u,d) == NeuronVals(u,respVali(u),2)
                    rectangle('Position', [posX(u,d)/(bin2)+round(preBase/bin2), (posY(u,d)*trialDivision-trialDivision)+0.5, window_size(2)/(bin2), trialDivision],...
                        'EdgeColor', 'r', 'LineWidth', 1,'LineStyle','-.');
                else
                    %                         rectangle('Position', [posX(u,d)/(bin2)+round(preBase/bin2), (posY(u,d)*trialDivision-trialDivision+1)/mergeTrials, window_size(2)/(bin2), trialDivision],...
                    %                         'EdgeColor', 'b', 'LineWidth', 0.5,'LineStyle','-.');
                end

            end


            xticklabels([])
            xlim([0 round(stimDur+preBase*2)/bin2])
            xticks([0 preBase/bin2:300/bin2:(stimDur+preBase)/bin2 (stimDur+preBase*2)/bin2])
            xticklabels([[0 preBase:300:(stimDur+preBase) round(stimDur+preBase*2)]]);


            %%%%%% Plot PSTH

            subplot(10,1,[7 8])

            MRhist = BuildBurstMatrix(goodU(:,u),round(p.t),round((directimesSorted(round(C(:,2)) == 0)-preBase)),round((stimDur+preBase*2)));


            [nT2,nN2,nB2]= size(MRhist);

            MRhist = squeeze(MRhist);

            spikeTimes = repmat([1:nB2],nT2,1);

            spikeTimes = spikeTimes(logical(MRhist));
            % Define bin edges (adjust for resolution)
            binWidth = 125; % 10 ms bins
            edges = [1:binWidth:round((stimDur+preBase*2))]; % Adjust time window as needed

            % Compute histogram
            psthCounts = histcounts(spikeTimes, edges);

            % Convert to firing rate (normalize by bin width)
            psthRate = (psthCounts / (binWidth * nT2))*1000;


            b=bar(edges(1:end-1), psthRate,'histc');
            b.FaceColor = 'k';
            b.FaceAlpha = 0.5;
            b.MarkerEdgeColor = "none";
            ylabel('Firing Rate (Hz)');
            xlim([0 round((stimDur+preBase*2)/100)*100])

            xticks([0 preBase:300:(stimDur+preBase*2) round((stimDur+preBase*2)/100)*100])
            xticklabels([])
            xline([preBase stimDur+preBase],'LineWidth',1.5)

            %set(gca, 'XLim', [-0.5, 1]); % Match raster plot x-axis


            %%%%PLot several trials one channel:

            chan = goodU(1,u);

            startTimes = directimesSorted(round(C(:,2)) == 0)-preBase;

            %startTimes = startTimes(end-9:end);

            window = round(stimDur+preBase*2);

            freq = "AP"; %or "LFP"

            type = "line"; %or heatmap

            subplot(10,1,[9 10])

            PlotRawData(fig,NP,chan,startTimes,window,freq,type,preBase);

            xticklabels([]) %%%plot example PV139-local neurons
            xlabel([])
            %xlim([0 round(stimDur+preBase*2)])


            xticklabels([-(preBase) 0:300:round((stimDur/100))*100 round((stimDur/100))*100 + preBase]./1000)


            %
            cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
            print(gcf, sprintf('%s-MovBall-%s-U%d-W%d-%dW-speed-500.pdf',NP.recordingName,orderNames{k},u,window_size(1),window_size(2)), '-dpdf', '-r300', '-vector');

            close


            %                 %%%Polar plot
            %                 cd(NP.recordingDir)
            %                 tuningCurve = (load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve)*1000; %convert to spikes/sec
            %                 theta = deg2rad(uDir); %linspace(0, 2*pi, size(tuningCurve,2)+1);  % 9 points for 8 bars (because it's circular)
            %                   % Remove the last value to avoid duplication of the first
            %                 pf = figure;
            %                 % Create the polar plot
            %                 polarplot([theta, theta(1)], [tuningCurve(eNeuron,:), tuningCurve(eNeuron,1)], '-o')
            %                 set(pf,"Color",'w')
            %                 ax = gca;
            %                 ax.ThetaTick = uDir;
            %                 title(sprintf('PolarPlot-U.%d-Unit-phy-%d',u,GoodU_or(u)));
            %                 cd(NP.recordingDir + "\Figs")
            %                 print(fig, sprintf('%s-MovBall-polarPlot-U%d.png',NP.recordingName,u),'-dpng');
            %                 close


        end
    end


end
% %%%%%%%%%%%%%%%%%%%%%%%% Position heatmap

for receptiveField=1 %%ReceptiveFieldMethod1 fixed delay

    if ReceptiveFieldFixedDelay ==1
        %%%

        %Get X and Y positions

        Xpos = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballTrajectoriesX'))));

        Ypos = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballTrajectoriesY'))));

        sizeN = length(unique(sizes));
        sizeX = size(Xpos);
        %%% X Y stucture = speed, offsets, directions, frames
        % A = [stimOn directions' offsets' sizes' speeds' orientations'] = Order of
        % categories (bigger divs to smaller divs.

        %%%Create a matrix with trials that have unique positions
        ChangePosX = zeros(sizeX(1)*sizeX(2)*sizeX(3)*sizeN*trialDivision,sizeX(4));
        ChangePosY = zeros(sizeX(1)*sizeX(2)*sizeX(3)*sizeN*trialDivision,sizeX(4));

        j=1;

        %For loop order determines the order of categories in ChangePosX = allTrials x nFrames
        for d = 1:sizeX(3) %directions
            for of = 1:sizeX(2) %offsets
                for sp = 1:sizeX(1) %speeds

                    ChangePosX(j:j+sizeN*trialDivision-1,:) = repmat(squeeze(Xpos(sp,of,d,:))',sizeN*trialDivision,1); %Size is not represented in X matrix.

                    ChangePosY(j:j+sizeN*trialDivision-1,:) = repmat(squeeze(Ypos(sp,of,d,:))',sizeN*trialDivision,1);

                    j = j+sizeN*trialDivision;

                end
            end
        end

        % Calculate the receptive field by multiplying normalized spike rate in frame by the image of each frame.

        goodNeurons = load(sprintf('pvalsBaselineBoot-1000-%s',NP.recordingName)).pvalsResponse;

        goodNeurons = find(goodNeurons<0.05);
        [Mb] = BuildBurstMatrix(goodU(:,goodNeurons),round(p.t),round((stims-preBase)),round(preBase)); %Baseline matrix plus

        Mb = mean(Mb,3); %mean across time bins

        spkRateBM = mean(Mb);

        %Original order
        % A = [stimOn directions' offsets' sizes' speeds' orientations']; %

        %Convolution order
        % A = [stimOn directions' offsets' speeds' orientations', sizes'];

        A = [A(:,1) A(:,2) A(:,3) A(:,5) A(:,6) A(:,2)]; % Reorganize matris
        [C indexS] = sortrows(A,[2 3 4 5 6]);

        %4. Sort directions:
        directimesSorted = C(:,1)';
        sizeV = C(:,6);

        delayResp = 200; %Based on unit 20 PV103_1
        msPerFarme= round(stimDur/sizeX(4));

        coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));

        [x, y] = meshgrid(1:coorRect(3)/10,1:coorRect(4)/10);

        sizesU = unique(sizeV);

        %Initialize matrices:
        matrixForNorm = zeros();

        RFu = zeros(sizeN,coorRect(4)/10,coorRect(3)/10,nN,"single");

        % Fill in the matrix with the counts of repeated coordinates
        SpeedSum = zeros(speedN,coorRect(4)/10,coorRect(3)/10,nN,"single");
        Stemp = zeros(coorRect(4)/10,coorRect(3)/10,nN,"single");
        DirecSum = zeros(direcN,coorRect(4)/10,coorRect(3)/10,nN,"single");
        Dtemp = Stemp;
        OrientSum = zeros(orientN,coorRect(4)/10,coorRect(3)/10,nN,"single");
        Otemp = Stemp;
        Si = 1;Di = 1;Oi = 1;

        for i = 1:trialDivision*sizeN:length(ChangePosX)

            %Matrix is sorted first by dir, then by size,
            [MtrialCat] = BuildBurstMatrix(goodU(:,goodNeurons),round(p.t),round((directimesSorted(i:i+trialDivision*sizeN-1))+delayResp),...
                round(mean(directimesSortedOff(i:i+trialDivision*sizeN-1)-directimesSorted(i:i+trialDivision*sizeN-1))));

            % Step 1: Reshape A to group elements to be averaged
            tempA = reshape(MtrialCat, [sizeN, size(MtrialCat,1)/sizeN, size(MtrialCat,2), size(MtrialCat,3)]);

            % Step 2: Calculate the mean along the 1st dimension
            tempB = squeeze(mean(tempA, 2));
            if length(size(tempB))<3
                tempB = reshape(tempB,[1,size(tempB,1),size(tempB,2)]);
            end

            % Step 3: Reshape the result to get the final matrix
            mSpkRate = tempB-spkRateBM; %T

            j = 1;
            tr = 1;
            while j <= length(MtrialCat) %size(ChangePosX,2)%length(MtrialCat) |

                matrix_nXnY = zeros(sizeN,coorRect(4)/10,coorRect(3)/10,"single");
                %matrixResp = zeros(coorRect(4),coorRect(3),nN,"single");

                centerX = ChangePosX(i,tr)/10;
                centerY = ChangePosY(i,tr)/10;
                radius = sizesU/2;

                % Calculate the distance of each point from the center
                distances = sqrt((x - centerX).^2 + (y - centerY).^2);

                for r =1:sizeN
                    % Set the values inside the circle to 1 (or any other value you prefer)
                    matrix_nXnY(r,distances <= radius(r)/10) =1;
                end

                frameSpkR = mean(mSpkRate(:,:,j:min(tr*msPerFarme,length(MtrialCat))),3);
                % Reshape frameSpkR to enable broadcasting
                ReshA = reshape(frameSpkR,[sizeN,1,1,nN]);

                %Repeat matrix_nXnY across the 4th dimension to match B's new shape
                ReshB = repmat(matrix_nXnY,[1,1,1,nN]);

                % Now, perform element-wise multiplication
                matrixResp = ReshB.*ReshA;

                if sum(isnan(matrixResp),'all')>1
                    1+1
                end

                %Feel individual params RF:
                %Order of organization = directions' offsets' speeds' orientations', sizes'

                if size(matrixResp,1) ~=1
                    matR = squeeze(mean(matrixResp));
                else
                    matR = squeeze(matrixResp);
                end
                if speedN > 1
                    Stemp = Stemp + matR;
                    if ismember((i-1)/(nT/direcN/offsetN/speedN), 1:speedN)
                        SpeedSum(Si,:,:,:) = squeeze(SpeedSum(Si,:,:,:))+Stemp;
                        Si = Si+1;
                        if Si == speedN
                            Si =1;
                            Stemp = zeros(coorRect(4)/10,coorRect(3)/10,nN,"single");
                        end
                    end
                end
                if direcN > 1
                    Dtemp = Dtemp + matR;
                    if ismember((i-1)/(nT/direcN),1:direcN) %nT/direcN:nT/direcN:nT
                        DirecSum(Di,:,:,:) = squeeze(DirecSum(Di,:,:,:))+Dtemp;
                        Di = Di+1;
                        if Di == direcN
                            Di =1;
                            Dtemp = zeros(coorRect(4)/10,coorRect(3)/10,nN,"single");
                        end
                    end
                end
                if orientN > 1
                    Otemp = Otemp + matR;
                    if ismember((i-1)/(nT/direcN/offsetN/speedN/orientN), 1:orientN)
                        OrientSum(Oi,:,:,:) = squeeze(OrientSum(Oi,:,:,:))+Otemp;
                        Oi = Oi+1;
                        if Oi == direcN
                            Oi =1;
                            Dtemp = zeros(coorRect(4)/10,coorRect(3)/10,nN,"single");
                        end
                    end
                end
                j = tr*msPerFarme;
                tr = tr+1;
                RFu = RFu+matrixResp;
                matrixForNorm = matrixForNorm+matrix_nXnY;
            end
        end
        normMatrix = repmat(matrixForNorm,[1,1,1,nN]).*reshape(spkRateBM,[1,1,1,length(spkRateBM)]);
        normRFu = RFu./normMatrix; %expected random rate

        cd(NP.recordingDir)
        save(sprintf('RFu_MovingBall-%s',NP.recordingName),'normRFu')
        save(sprintf('NormMatrix_MovingBall-%s',NP.recordingName),'normMatrix')
        save(sprintf('RFuSpeed_MovingBall-%s',NP.recordingName),'SpeedSum')
        save(sprintf('RFuOrient_MovingBall-%s',NP.recordingName),'OrientSum')
        save(sprintf('RFuDirec_MovingBall-%s',NP.recordingName),'DirecSum')


    end
end

for depthPlots = 1
    if depthPlot == 1

        % Depth plot
        if ~exist('FigSig','var')

            FigSig = figure;
        end


        pvals = load(sprintf('pvalsBaselineBoot-%d-%s.mat',N_bootstrap,NP.recordingName)).pvalsResponse;

        goodUR = goodU(:,pvals<0.005);

        respSpR = max(NeuronVals(pvals<0.005,:,1),[],2)*(1000/duration); %spikes per duration in window ms converted to spikes per second

        % Normalize responses to range from 0 (black) to 1 (white)
        colors = repmat((respSpR - min(respSpR)) / (max(respSpR) - min(respSpR)),[1,3]);

        [Coor] = neuronLocation(NP,data(ex,:),goodUR,0,0,0);

        figure
        scatter(respSpR,Coor{1,3}(2,:))
        set(gca, 'YDir', 'reverse');



        NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName));
        tuningCurve = load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve;
        %
        %     Theta = load(sprintf('Angle-prefer-%s',NP.recordingName)).Theta;
        %     TI= load(sprintf('Tuning-Index-%s',NP.recordingName)).L;

        Lo



        respU = string(data.ResponseU(ex));


        % strlength(respU);
        if strlength(respU)>0

            strSplit = strsplit(respU, ','); % split the string
            ResponsiveN = cellfun(@str2num, strSplit); % convert to numbers

            ResponsiveDepth =  verticalDepthU(ResponsiveN);

            [CM,CMP,h, cmapR_theta]=colormap2D('plotPolarExample',1); %create colormap

            [nRC,~,nTC]=size(cmapR_theta);

            rhoR = TI(ResponsiveN);

            thetaR = Theta(ResponsiveN)+90; %sum 90 degrees because 0 is north (90)

            rhoR(rhoR == 0) = 0.001;

            thetaR(thetaR == 0) = 0.001;

            divFactor = max(unique(thetaR));


            hold on;
            for i =1:length(thetaR)
                colorsN = cmapR_theta( ceil((thetaR(i)/(divFactor))*nRC), : , ceil(rhoR(i)/nTC));
                scatter(1+(ex-x)*5, ResponsiveDepth(i),'filled','CData',colorsN,'MarkerEdgeColor','k', 'LineWidth',0.5, 'SizeData',40)

            end

            ylabel('Unit depth (um)')
            xlabel('Insertions')
        else
            x = x+1;
        end



    end
end

% %%% Convolution
for convNeuron = 1
%     if ismember(ex,GoodRecordingsRF)
%         ReceptiveFieldConvolutions =1;
%         repeatConv =1;    
%     else
%         ReceptiveFieldConvolutions =0;
%         repeatConv =0; 
%     end

    if ReceptiveFieldConvolutions ==1
        cd(NP.recordingDir)

        pvalTi= load(sprintf('pvalsBaselineBoot-1000-%s',NP.recordingName)).pvalsResponse;

        respU = find(pvalTi <0.05);

        if ~isempty(respU)


        if ~isfile(sprintf('RFuC-%s.mat',NP.recordingName)) || repeatConv%exist
            %%A. CONVOLUTION SETUP
            %%%%%%%%%%Load responsive neurons

            %%%%%%%%%%Load X and Y ball positions
            nFrames = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nFrames'))));
            framesFast = min(unique(nFrames)); %Select number of frames from faster speed

            Xpos = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballTrajectoriesX'))));
            Ypos = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballTrajectoriesY'))));
           
            
            if size(Xpos,1) >1
                Xpos = Xpos(2,:,:,1:framesFast); %Select number of frames from faster speed
                Ypos = Ypos(2,:,:,1:framesFast);
                %%2 speeds
            end

            

            sizeN = length(unique(sizes));
            sizeX = size(Xpos());
            %%% X Y stucture = speed, offsets, directions, frames
            % A = [stimOn directions' offsets' sizes' speeds' orientations'] = Order of
            % categories (bigger divs to smaller divs.

            %%%Create a matrix with trials that have unique positions
            ChangePosX = zeros(sizeX(1)*sizeX(2)*sizeX(3)*sizeN*trialDivision,sizeX(4));
            ChangePosY = zeros(sizeX(1)*sizeX(2)*sizeX(3)*sizeN*trialDivision,sizeX(4));
            j=1;
            sizeX = size(Xpos);

            %For loop order determines the order of categories in ChangePosX = allTrials x nFrames
            for d = 1:sizeX(3) %directions
                for of = 1:sizeX(2) %offsets
                    for sp = 1:sizeX(1) %speeds

                        ChangePosX(j:j+sizeN*trialDivision-1,:) = repmat(squeeze(Xpos(sp,of,d,:))',sizeN*trialDivision,1); %Size is not represented in X matrix.

                        ChangePosY(j:j+sizeN*trialDivision-1,:) = repmat(squeeze(Ypos(sp,of,d,:))',sizeN*trialDivision,1);

                        j = j+sizeN*trialDivision;

                    end
                end
            end

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if noEyeMoves %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                  % EyePositionAnalysis
                  % Create spike Sums with NaNs when the eye is not present.
                  
        
            file = dir (NP.recordingDir);
            filenames = {file.name};
            files= filenames(contains(filenames,"timeSnipsNoMov"));
            cd(NP.recordingDir)
            %Run eyePosition Analysis to find no movement timeSnips
            if XYdivision %%%Add time snips for x division and for y division (for loop 1:4)
                xfiles = filenames(contains(filenames,"timeSnipsNoMov")& contains(filenames,"X-change"));
                timeSnipsX = load(xfiles{numel(xfiles)}).timeSnips;
                Times{1} = timeSnipsX;
                
                yfiles = filenames(contains(filenames,"timeSnipsNoMov")& contains(filenames,"Y-change"));
                timeSnipsY = load(yfiles{numel(yfiles)}).timeSnips;
                Times{2} = timeSnipsY;
                

            else

                timeSnips = load(files{numel(files)}).timeSnips;
                timeSnipsMode = timeSnips(:,timeSnips(3,:) == mode(timeSnips(3,:)));
                Times{1} = timeSnips;
            end
% 
%             A = [stimOn directions' offsets' sizes' speeds' orientations'];
%             [C indexS] = sortrows(A,[2 3 4 5 6]);
% 
%             %find(-stimOn+stimOff>3000)
% 
% 
%             %4. Sort directions:
%             directimesSorted = C(:,1)';

%            %Original order
           % A = [stimOn directions' offsets' sizes' speeds' orientations']; %

           %When dealing with different speeds, save different stim durs, and create M for each speed

           %[C indexS] = sortrows(A,[2 3 4 5 6]);

           %Convolution order
            %A = [stimOn directions' offsets' speeds' orientations', sizes'];

           A = [A(:,1) A(:,2) A(:,3) A(:,5) A(:,6) A(:,4)]; % Reorganize matrix
           [C indexS] = sortrows(A,[2 3 4 5 6]);

           %4. Sort directions:
           directimesSorted = C(:,1)';
           sizeV = C(:,6);

           %%%%% Asign trial index accordying to pressence of trial within
           %%%%% eye timesnip:

           if ModeQuadrant %%I can make times for each bin, if the time of each bin is within the quadrant, then assign.
               selecInd = [];
               for i = 1:size(timeSnipsMode,2)

                   %Find stimOns and offs that are between each timeSnip
                   selecInd = [selecInd find(directimesSorted>=timeSnipsMode(1,i) & directimesSorted<(timeSnipsMode(2,i)-stimDur))];
               end

               IndexQ{1} = selecInd;
               IndexDiv{1} = IndexQ;

           else %%%%%%% get indexes of directed times sorted for

               for t = 1:numel(Times)

                   timeSnips = Times{t};
                   for j = unique(timeSnips(3,:))
                       selecInd = [];
                       timeSnipsQ = timeSnips(:,timeSnips(3,:) ==j);
                       for i = 1:size(timeSnipsQ,2)

                           %Find stimOns and offs that are between each timeSnip
                           selecInd = [selecInd find(directimesSorted>=timeSnipsQ(1,i) & directimesSorted<(timeSnipsQ(2,i)-stimDur))];

                       end

                       IndexQ{j} = selecInd;
                   end
                   IndexDiv{t} = IndexQ;
               end

           end

           coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));
           reduceFactor = min([20 min(sizeV)]); %has to be bigger than the smallest ball size
           redCoorX = round(coorRect(3)/reduceFactor);
           redCoorY = round(coorRect(4)/reduceFactor);

           [x, y] = meshgrid(1:redCoorX,1:redCoorY);

           sizesU = unique(sizeV);

           x = fliplr(x);

            %%%%%% Add spikes accordying to number of frames.

            msPerFarme= stimDur/sizeX(4);

            duration = 300;
            
            if useZscore %%%%%Select correct timestamps.
                Mr = load(sprintf('ZscoreRaster-%d-%s',N_bootstrap,NP.recordingName)).ZscoreRaster;
            else
                [Mr] = BuildBurstMatrix(goodU,round(p.t/msPerFarme),round((directimesSorted)/msPerFarme),round((stimDur)/msPerFarme));
                [Mbase] = BuildBurstMatrix(goodU,round(p.t/msPerFarme),round((directimesSorted-duration)/msPerFarme),round((duration)/msPerFarme));
            end

            [trials neurons bins] = size(Mr());
            
            spikeSumsDiv = cell(1,numel(IndexDiv));
            spikeSumsQ = cell(1,numel(IndexQ));

            for t = 1:numel(IndexDiv) %%%Divison, horizontal or vertical

                for i = 1:numel(IndexQ) %%%% quadrant, e.g. right, left

                    spikeSum = zeros(size(Mr,1),length(respU),sizeX(4),'single');

                    %spikeSum(IndexDiv{t}{i},:,:) = Mr(IndexDiv{t}{i},respU,:) - mean(Mbase(IndexDiv{t}{i},respU,:),3); %Substract baseline per trial

                    spikeSum(IndexDiv{t}{i},:,:) = Mr(IndexDiv{t}{i},respU,:); %Substract baseline per trial

                    spikeSum = (spikeSum./msPerFarme).*1000; %Convert to spikes per second
                    
                    %Previous method:
%                    spikeSum = zeros(size(Mr,1),length(respU),sizeX(4),'single');
%
%                     for u = 1:length(respU)
% 
%                         Mu = squeeze(Mr(:,respU(u),:));
% 
%                         j= 1;
% 
%                         for f = 1:sizeX(4)
% 
%                             spikeSum(IndexDiv{t}{i},u,f) = mean(Mu(IndexDiv{t}{i},1*j:min(f*msPerFarme,length(Mu))),2);
% 
%                             j = f*msPerFarme;
%                         end
%                     end

                    spikeSumsQ{i} = spikeSum;
                end

                spikeSumsDiv{t} = spikeSumsQ;
            end

            cd(NP.recordingDir)


            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            else %%%%%%%%%%% Eye moves: %%%%%%%%%%%%%%%%%%%%%%

            
            %Original order
           % A = [stimOn directions' offsets' sizes' speeds' orientations']; %

            %Convolution order
           % A = [stimOn directions' offsets' speeds' orientations', sizes'];

            A = [A(:,1) A(:,2) A(:,3) A(:,5) A(:,6) A(:,4)]; % Reorganize matrix
            [C indexS] = sortrows(A,[2 3 4 5 6]);

            %4. Sort directions:
            directimesSorted = C(:,1)';
            sizeV = C(:,6);

            coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));

            reduceFactor = min([20 min(sizeV)]); %has to be bigger than the smallest ball size

            redCoorX = round(coorRect(3)/reduceFactor);
            redCoorY = round(coorRect(4)/reduceFactor);

            [x, y] = meshgrid(1:redCoorX,1:redCoorY);

            sizesU = unique(sizeV);

            x = fliplr(x);

            msPerFarme= stimDur/sizeX(4);
            
            if useZscore %%%%%Select correct timestamps.
                ZscoreRaster = load(sprintf('ZscoreRaster-%d-%s',N_bootstrap,NP.recordingName)).ZscoreRaster;
            else
                [Mr] = BuildBurstMatrix(goodU(:,respU),round(p.t/msPerFarme),round((directimesSorted)/msPerFarme),round((stimDur)/msPerFarme));
                [Mbase] = BuildBurstMatrix(goodU(:,respU),round(p.t/msPerFarme),round((directimesSorted-preBase)/msPerFarme),round((preBase)/msPerFarme));
            end

            [trials neurons bins] = size(Mr);
% 
%             for u = 1:length(respU)
% 
%                 Mu = squeeze(Mr(:,respU(u),:));
% 
%                 j= 1;
% 
%                 for f = 1:sizeX(4)
% 
%                     spikeSum(:,u,f) = mean(Mu(:,1*j:min(f*msPerFarme,length(Mu))),2);
% 
%                     j = f*msPerFarme;
%                 end
%             end
%             cd(NP.recordingDir)
%            boot_means = load(sprintf('MovBall-Base-Boot-1000-%s',NP.recordingName)).boot_means;

            %%%Zscore of spikes per frame. Convert everything to spikes per
            %%%second beforehand. 
%             meanBase = mean(Mbase(:,respU,:),[1 3]);
%             N_bootstrap = 1000;
%             nN= length(respU);
%             substractor = reshape(mean(boot_means(:,respU)*(1000/duration)),[1, nN, 1]); 
%             denominator = reshape(std(boot_means(:,respU)*(1000/duration))+1/(N_bootstrap*trialDivision),[1,nN,1]);
%             spikeSum = (spikeSum.*(1000/msPerFarme)-substractor)./denominator; %convert into spike rate.
           
            spikeSum = (Mr./msPerFarme).*1000; %Convert to spikes per second


%             ZscoreRaster = load(sprintf('ZscoreRaster-%d-%s',N_bootstrap,NP.recordingName)).ZscoreRaster;
% 
%             %bin Z-score accordying to frame number:
% 
% 
%             for u = 1:length(respU)
% 
%                 Mu = squeeze(ZscoreRaster(:,respU(u),:));
% 
%                 j= 1;
% 
%                 for f = 1:sizeX(4)
% 
%                     spikeSum(:,u,f) = mean(Mu(:,1*j:min(f*msPerFarme,length(Mu))),2);
% 
%                     j = f*msPerFarme;
%                 end
%             end

            spikeSumsDiv{1}{1} = spikeSum;

            IndexDiv = 1;

            IndexQ = 1;

            end %%%%%%%%%%%%%%%%%%%%%%%%%%%%% end eyeMovement initialization
%
          
%             
%             %%%Test convolution
%             spikeSumArt = zeros(size(spikeSum));
%             spikeSumArt([1:10 91:100 181:190 271:280],7,end-20:end)=1;%Spiking at the end of first offset
%             %spikeSumArt(nT/4+1:nT/4+15,7,end-20:end) =1;%
%             spikeSum = spikeSumArt;
%             figure;imagesc(squeeze(spikeSum(:,7,:)));colormap(flipud(gray(64)));
%             yline(trialDivision*sizeN:trialDivision*sizeN:size(spikeSumArt,1));
%             yline(trialDivision*offsetN*sizeN:trialDivision*offsetN*sizeN:size(spikeSumArt,1)-1,'r','LineWidth',5)

            %%A. CONVOLUTION. Runs the spike train across the stimulus videos, to
            %%extract noisy receptive field (because spike responses are noise).
            if noEyeMoves
                Q = unique(timeSnips(3,:));
            else
                Q = 1;
            end

            %%%%%%%%%%%%%%%%%%%% Shuffle raster before convolution in order
            %%%%%%%%%%%%%%%%%%%% to calculate tuning index
            %%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%

            nShuffle =50;

            Raster = BuildBurstMatrix(goodU(:,respU),round(p.t/msPerFarme),round((directimesSorted)/msPerFarme),round((stimDur)/msPerFarme));
            
            shuffledData = zeros(size(Raster,1), size(Raster,2), size(Raster,3), nShuffle);

            for i =1:nShuffle

                % Shuffle along the first dimension
                idx1 = randperm(size(Raster,1));

                % Shuffle along the third dimension
                idx3 = randperm(size(Raster,3));

                shuffledData(:,:,:,i) = Raster(idx1, :, idx3);

            end

            tic


         for t = 1:numel(IndexDiv)

             for q = 1:numel(IndexQ)

                 %%%%Initialize 5D matrices
                 RFuShuff = zeros(redCoorY,redCoorY,sizeX(4),length(respU),"single");
                 RFu = zeros(redCoorY,redCoorY,sizeX(4),length(respU),"single");

                 RFuSpeed = zeros(speedN,redCoorY,redCoorY,sizeX(4),length(respU),"single");
                 RFuDir = zeros(direcN,redCoorY,redCoorY,sizeX(4),length(respU),"single");
                 RFuSize = zeros(sizeN,redCoorY,redCoorY,sizeX(4),length(respU),"single");
                 RFuDirSize = zeros(direcN,sizeN,redCoorY,redCoorY,sizeX(4),length(respU),"single");


                 Uspeed = unique(speeds);
                 Udir = unique(directions);
                 Usize = unique(sizes);
                 NormVideo = zeros(redCoorY,redCoorY,sizeX(4),'single');

                 spikeSum = spikeSumsDiv{t}{q};
                 for i = 1:trialDivision:trials

                     videoTrials = zeros(redCoorY,redCoorY,sizeX(4),'single');
                     for j = 1:sizeX(4) %%Calculate video of unique trials


                         xyScreen = zeros(redCoorY,redCoorX,"single");
                         %matrixResp = zeros(coorRect(4),coorRect(3),nN,"single");

                         centerX = ChangePosX(i,j)/reduceFactor;
                         centerY = ChangePosY(i,j)/reduceFactor;
                         radius = sizeV(i)/2;

                         % Calculate the distance of each point from the center
                         distances = sqrt((x - centerX).^2 + (y - centerY).^2);

                         % Set the values inside the circle to 1 (or any other value you prefer)
                         xyScreen(distances <= radius/reduceFactor+0.5) = 1;

                         videoTrials(:,:,j) = xyScreen(:,(redCoorX-redCoorY)/2+1:(redCoorX-redCoorY)/2+redCoorY);
                         %figure;imagesc(xyScreen);
                     end

                     %spikeMean = mean(spikeSum(i:i+trialDivision-1,:,:)-spkRateBM);
                     %Normalize spike mean by number of spikes.
                     spikeMean = mean(spikeSum(i:i+trialDivision-1,:,:),'omitnan');
                     Co = zeros(redCoorY,redCoorY,sizeX(4),length(respU),'single');
                     CoShuff = zeros(redCoorY,redCoorY,sizeX(4),length(respU),nShuffle,'single');

                     spikeMeanShuff = mean(shuffledData(i:i+trialDivision-1,:,:,:),'omitnan');


                     for u = 1:length(respU)  

                         if C(i,2) ==uDir(1) || C(i,2) ==uDir(3) %%%%%Rotate convolution if angle is up or down in order to reflect true video. 
                            Co(:,:,:,u) =rot90(convn(videoTrials,spikeMean(:,u,:),'same'),2);
                         else
                             Co(:,:,:,u) =convn(videoTrials,spikeMean(:,u,:),'same');
                         end

                         %%%Convolve the shuffled rasters

                         for s = 1:nShuffle

                              if C(i,2) ==uDir(1) || C(i,2) ==uDir(3)
                                  CoShuff(:,:,:,u,s) = rot90(convn(videoTrials,spikeMeanShuff(:,u,:,s),'same'),2);
                              else
                                  CoShuff(:,:,:,u,s) = convn(videoTrials,spikeMeanShuff(:,u,:,s),'same');
                              end

                         end                        

                     end

                   
                     RFuDir(Udir == C(i,2),:,:,:,:) = squeeze(RFuDir(Udir == C(i,2),:,:,:,:))+Co./(nT/direcN/trialDivision);
                     RFuDirSize(Udir == C(i,2),Usize == C(i,6),:,:,:,:) = squeeze(RFuDirSize(Udir == C(i,2),Usize == C(i,6),:,:,:,:))+Co./(nT/direcN/trialDivision);
                     RFuSize(Usize == C(i,6),:,:,:,:) = squeeze(RFuSize(Usize == C(i,6),:,:,:,:))+Co./(nT/sizeN/trialDivision);
                     RFuSpeed(Uspeed == C(i,4),:,:,:,:) = squeeze(RFuSpeed(Uspeed == C(i,4),:,:,:,:))+Co./(nT/speedN/trialDivision);

                     NormVideo = NormVideo+videoTrials;%.*reshape(spkRateBM,[1 1 1 length(spkRateBM)]);

                     RFu = RFu+Co./(nT/trialDivision);
                     RFuShuff = RFuShuff+CoShuff./(nT/trialDivision);
                 end

                 toc
                 %implay(squeeze(RFu(:,:,:,1)));
                 %implay(videoTrials)

                 %%%%%%%%%% Normalization parameters
                 L = size(spikeSum,3);
                 time_zero_index = ceil(L / 2);


                 nN =length(respU);

                 %%%%%% Normalize and select delay
                 delay = 250;

                 %             RFuNorm = RFu;%squeeze(RFuSize(ceil(sizeN/2),:,:,:,:))./reshape(normMatrixMean,[size(normMatrixMean,1) size(normMatrixMean,2) 1 nN]);
                 %             RFuST = squeeze(RFuNorm(:,:,time_zero_index+round(delay/msPerFarme),:));
                 %             RFuNormDir = (RFuDir-normMatrixMean)./normMatrixSTD;
                 %             %%%%%RFuSTDir = squeeze(RFuNormDir(:,:,:,time_zero_index+round(delay/msPerFarme),:));
                 RFuSTDir =  squeeze(RFuDir(:,:,:,time_zero_index+round(delay/msPerFarme),:));

                 RFuSTDirSize =  reshape(RFuDirSize(:,:,:,:,time_zero_index+round(delay/msPerFarme),:),... %Reshape to eliminate frame component.
                     [size(RFuDirSize,1),size(RFuDirSize,2),size(RFuDirSize,3),size(RFuDirSize,4),size(RFuDirSize,6)]);

                 RFuST = squeeze(RFu(:,:,time_zero_index+round(delay/msPerFarme),:));

                 RFuShuffST = squeeze(RFuShuff(:,:,time_zero_index+round(delay/msPerFarme),:,:));                
               

                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Apply 2d Gaussian to receptive fields
                 %TwoDGaussian = fspecial('gaussian',floor(size(RFuSTDirSize,3)/9),5);

                 TwoDGaussian = fspecial('gaussian',floor(size(RFuSTDirSize,3)/(9/2)),redCoorY/offsetN); %increase size of gaussian by 100%.
                 %RFuSTDirFilt = zeros(size(RFuSTDir));
                 
                 RFuSTDirSizeFilt = zeros(size(RFuSTDirSize));

               
                 for d = 1:size(RFuSTDir,1) %dirs
                     for s = 1:size(RFuSTDirSize,2) %Size
                         parfor ui =1:size(RFuSTDir,4) %units

                             slice = squeeze(RFuSTDirSize(d,s,:,:,ui));

                             slicek = conv2(slice,TwoDGaussian,'same');

                             RFuSTDirSizeFilt(d,s,:,:,ui) =slicek;

                         end
                     end

                 end


                 names = {'X','Y'};
                 if noEyeMoves
                     save(sprintf('NEM-RFuSTDirSizeFilt-Q%d-Div-%s-%s',q,names{t},NP.recordingName),'RFuSTDirSizeFilt','-v7.3')
                     save(sprintf('NEM-RFuSTDirSize-Q%d-Div-%s-%s',q,names{t},NP.recordingName),'RFuSTDirSize','-v7.3')
                     save(sprintf('NEM-RFuST-Q%d-Div-%s-%s',q,names{t},NP.recordingName),'RFuST','-v7.3')
                     save(sprintf('NEM-RFuShuffST-Q%d-Div-%s-%s',q,names{t},NP.recordingName),'RFuShuffST','-v7.3')

                     %save(sprintf('NEM-spikeSums-%s',NP.recordingName),'spikeSums','-v7.3')
                     %                 save(sprintf('NEM-RFuSelecTime-%s',NP.recordingName),'RFuST','-v7.3')
                     %                 save(sprintf('NEM-RFuSelecTimeMask-%s',NP.recordingName),'RFuST','-v7.3')
                     %                 save(sprintf('NEM-RFuSelecTimeD-%s',NP.recordingName),'RFuSTDir','-v7.3')
                     %                 save(sprintf('NEM-RFuSTDirFilt-%s',NP.recordingName),'RFuSTDirFilt','-v7.3')
                     %                 save(sprintf('NEM-RFuSelecTimeMask-%s',NP.recordingName),'RFuSTmaskD','-v7.3')
                     %                 save(sprintf('NEM-RFuC-%s',NP.recordingName),'RFu','-v7.3')
                     %                 save(sprintf('NEM-RFuDirC-%s',NP.recordingName),'RFuDir','-v7.3')
                     %                 save(sprintf('NEM-RFuSizeC-%s',NP.recordingName),'RFuSize','-v7.3')
                     %                 save(sprintf('NEM-RFuSpeedC-%s',NP.recordingName),'RFuSpeed','-v7.3')
                     %                 save(sprintf('NormVideo-%s',NP.recordingName),'NormVideo','-v7.3')

                 else
                     save(sprintf('RFuSTDirSizeFilt-%s',NP.recordingName),'RFuSTDirSizeFilt','-v7.3')
                     save(sprintf('RFuSTDirSize-Q%d-Div-%s-%s',q,names{t},NP.recordingName),'RFuSTDirSize','-v7.3')
                     save(sprintf('RFuST-Q%d-Div-%s-%s',q,names{t},NP.recordingName),'RFuST','-v7.3')
                     save(sprintf('RFuShuffST-Q%d-Div-%s-%s',q,names{t},NP.recordingName),'RFuShuffST','-v7.3')
%                      save(sprintf('RFuNorm-%s',NP.recordingName),'RFuNorm','-v7.3')
%                      save(sprintf('RFuSelecTime-%s',NP.recordingName),'RFuST','-v7.3') %%Not normalized
%                      save(sprintf('RFuSelecTimeMask-%s',NP.recordingName),'RFuST','-v7.3')
%                      save(sprintf('RFuSelecTimeD-%s',NP.recordingName),'RFuSTDir','-v7.3') %%Normalized
%                      save(sprintf('RFuSTDirFilt-%s',NP.recordingName),'RFuSTDirFilt','-v7.3')
%                      %save(sprintf('RFuSelecTimeMask-%s',NP.recordingName),'RFuSTmaskD','-v7.3')
%                      save(sprintf('RFuC-%s',NP.recordingName),'RFu','-v7.3')
%                      save(sprintf('RFuDirC-%s',NP.recordingName),'RFuDir','-v7.3')
%                      save(sprintf('RFuSizeC-%s',NP.recordingName),'RFuSize','-v7.3')
%                      save(sprintf('RFuSpeedC-%s',NP.recordingName),'RFuSpeed','-v7.3')
%                      save(sprintf('NormVideo-%s',NP.recordingName),'NormVideo','-v7.3')

                 end
             end
         end


% %Test with figures:
%             figure;imagesc(squeeze(RFuSTmaskD(d,:,:,2)));
% 
% 
% % % 
% % %             testRFU = squeeze(RFuSize(1,:,:,:,4));
% % % 
% % %             max(testRFU,[],'all')
% % 
% % % 
%            
% % % 
% %              implay(RFuST(:,:,:,1));
% %              implay(NormVideo)

        else
            %RFuSTmask = load(sprintf('RFuSelecTime-%s',NP.recordingName)).RFuSTmask;
            %             RFuNorm = load(sprintf('RFuNorm-%s',NP.recordingName)).RFuNorm;
            %             RFu = load(sprintf('RFuC-%s',NP.recordingName),'RFu').RFu;
            %             RFuDir = load(sprintf('RFuDirC-%s',NP.recordingName),'RFuDir').RFuDir;
            %             RFuSize =  load(sprintf('RFuSizeC-%s',NP.recordingName),'RFuSize').RFuSize;
            %             RFuSpeed = load(sprintf('RFuSpeedC-%s',NP.recordingName),'RFuSpeed').RFuSpeed;
            if noEyeMoves
                RFuSTDir = load(sprintf('NEM-RFuSelecTimeD-%s',NP.recordingName)).RFuSTDir;
            else
                RFuSTDir = load(sprintf('RFuSelecTimeD-%s',NP.recordingName)).RFuSTDir;
            end
            
            %figure;imagesc(RFuSTmask(:,:,2));

        end
 
        if plotRF

%             eNeuron = find(ismember(respU,[29]));%respU(respU == selecN{1}(2,selecN{1}(1,:)==ex));
%             %figure;imagesc(squeeze(RFu(:,:,88,2)))
            eNeuron =40;
            %Parameters
            eye_to_monitor_distance = 21.5; % Distance from eye to monitor in cm
            pixel_size = 33/(1080/reduceFactor); % Size of one pixel in cm (e.g., 25 micrometers)
            monitor_resolution = [redCoorX, redCoorY]; % Width and height in pixels
            [theta_x theta_y] = pixels2eyeDegrees(eye_to_monitor_distance,pixel_size,monitor_resolution);

            for u = eNeuron
                fig = tiledlayout(direcN/2,direcN/2);
                
                for d = 1:direcN
                    if d ==1 || d==3
                    nexttile;imagesc(rot90(squeeze(RFuSTDir(d,:,:,u)),2));c = colorbar;caxis([0 max(RFuSTDir(:,:,:,u),[],'all')]);
                    else
                    nexttile;imagesc((squeeze(RFuSTDir(d,:,:,u))));c = colorbar;caxis([0 max(RFuSTDir(:,:,:,u),[],'all')]);    
                    end
                    colormap('jet')
                    title(string(uDir(d)))
                    title(c,'Z-score')
                    xlim([(redCoorX-redCoorY)/2 (redCoorX-redCoorY)/2+redCoorY])
                    
                    xt = xticks;%(linspace((redCoorX-redCoorY)/2,(redCoorX-redCoorY)/2+redCoorY,offsetN*2));
                    xticklabels(round(theta_x(1,xt)))
                    yt = yticks;
                    yticklabels(round(theta_y(yt,1)))
                    xlabel('X degrees')
                    ylabel('Y degrees')
                    
                end

                set(gcf,'Color','w')
                fig.Position = [0.13 0.114444438994877 0.722095239372481 0.810555561005123];

               
                colorbarlims = [0 max(RFuSTDir(:,:,:,eNeuron),[],'all')];
                cd(NP.recordingDir)
                save(sprintf('%s-Unit-%d-MovBall-RFlims-Dirs',NP.recordingName,respU(u)),'colorbarlims')
                cd(NP.recordingDir + "\Figs")
                if noEyeMoves
                print(gcf,sprintf('%s-NEM-Unit-%d-MovBall-RF-Dirs.png',NP.recordingName,respU(u)),'-dpng')
                else
                print(gcf,sprintf('%s-Unit-%d-MovBall-RF-Dirs.png',NP.recordingName,respU(u)),'-dpng')
                end

            end

        end
%
        if calculateEntropy

            entropies = zeros(1,length(respU));

            for u = 1:length(respU)


                M = squeeze(RFuSTDir(:,:,:,u));

                % Find the maximum value and its index
                [maxValue, linearIndex] = max(M(:));

                % Convert the linear index to subscripts to find the slice
                [sliceIndex, ~, ~] = ind2sub(size(M), linearIndex);

                % Extract the corresponding 10x10 slice
                Mred = squeeze(M(sliceIndex, :, (size(M,3)-size(M, 2))/2+1: size(M,3)-(size(M,3)-size(M, 2))/2));

                %%Prepare matrix into the same style as rectGrid

                % Original 50x100 matrix
                %M = rand(50, 100); % Replace with your matrix

                % Number of blocks (offsetN)
                offsetN = 9;

                % Define edges for rows and columns
                rowEdges = round(linspace(1, size(Mred, 1) + 1, offsetN + 1));
                colEdges = round(linspace(1, size(Mred, 2) + 1, offsetN + 1));

                % Initialize the resulting matrix
                reducedMatrix = zeros(offsetN, offsetN);

                % Compute mean for each block, omitting NaNs
                for i = 1:offsetN
                    for j = 1:offsetN
                        % Extract block using the calculated edges
                        block = Mred(rowEdges(i):rowEdges(i+1)-1, colEdges(j):colEdges(j+1)-1);
                        % Compute mean, omitting NaNs
                        reducedMatrix(i, j) = mean(block(:), 'omitnan');
                    end
                end

                % Normalize to create a probability distribution
                reducedMatrix = reducedMatrix / sum(reducedMatrix(:));

                % Convert to an image-like format and calculate entropy
                % (scale to [0, 1] for compatibility with `entropy`)
                P_scaled = mat2gray(reducedMatrix);
                figure;imagesc((P_scaled));
                set(gca, 'YDir', 'reverse', 'XDir', 'reverse');
                colorbar;

                entropies(u) = entropy(P_scaled);
            end


            cd(NP.recordingDir)
            sign = '0.05';

            if noEyeMoves

                save(sprintf('NEM-Entropies-MB-RF-respU-%s-%s.mat',sign,NP.recordingName),'entropies')

            else

                save(sprintf('Entropies-MB-RF-respU-%s-%s.mat',sign,NP.recordingName),'entropies')

            end
        

        end

        end

        end
end


% %%%%%%%%%%%%Spatial tuning
for spatun =1
    if spatialTuning ==1

        cd(NP.recordingDir)

        pvalTi= load(sprintf('pvalsBaselineBoot-1000-%s',NP.recordingName)).pvalsResponse;

        respU = find(pvalTi <0.005);

        %%%Divide screen into rect grid coordinates and calculate the
        %%%response of moving ball for each of this coordinates
        for extractRG =1
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
        offsetsR = [];
        sizesR = [];
        seqMatrix = [];


        if ~contains(data.VS_ordered(ex),"RG")
            %disp()
            w= sprintf('No rectangle grid files where found in %s. Skipping into next experiment.',NP.recordingName);
            warning(w)
            %continue
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
                    sizesR = [sizesR cell2mat(rect.VSMetaData.allPropVal(find(strcmp(rect.VSMetaData.allPropName,'tilingRatios'))))];
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
            sizesR = repmat(cell2mat(rect.VSMetaData.allPropVal(5)),1,length(seqMatrix));
        end

        screenSideO = rect.VSMetaData.allPropVal{find(strcmp(rect.VSMetaData.allPropName,'rect'))};
        screenSide = screenSideO(4);
        gridSize = length(unique(positionsMatrix));
        squareSize = screenSide / gridSize; % Size of each square
        end


        %%%%% Spikes and moving ball positions %%%%%%%
        for movBallpos =1

            %%%%%%%%%%Load X and Y ball positions
            Xpos = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballTrajectoriesX'))));

            Ypos = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballTrajectoriesY'))));

            sizeN = length(unique(sizes));
            sizeX = size(Xpos);
            %%% X Y stucture = speed, offsets, directions, frames
            % A = [stimOn directions' offsets' sizes' speeds' orientations'] = Order of
            % categories (bigger divs to smaller divs.

            %%%Create a matrix with trials that have unique positions
            ChangePosX = zeros(sizeX(1)*sizeX(2)*sizeX(3)*sizeN*trialDivision,sizeX(4));
            ChangePosY = zeros(sizeX(1)*sizeX(2)*sizeX(3)*sizeN*trialDivision,sizeX(4));
            j=1;
            sizeX = size(Xpos);

            %For loop order determines the order of categories in ChangePosX = allTrials x nFrames
            for d = 1:sizeX(3) %directions
                for of = 1:sizeX(2) %offsets
                    for sp = 1:sizeX(1) %speeds

                        ChangePosX(j:j+sizeN*trialDivision-1,:) = repmat(squeeze(Xpos(sp,of,d,:))',sizeN*trialDivision,1); %Size is not represented in X matrix.

                        ChangePosY(j:j+sizeN*trialDivision-1,:) = repmat(squeeze(Ypos(sp,of,d,:))',sizeN*trialDivision,1);

                        j = j+sizeN*trialDivision;

                    end
                end
            end

            A = [stimOn directions' offsets' speeds' orientations', sizes'];
            [C indexS] = sortrows(A,[2 3 4 5 6]);

            %4. Sort directions:
            directimesSorted = C(:,1)';
            sizeV = C(:,6);

            coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));

            reduceFactor = min([20 min(sizeV)]); %has to be bigger than the smallest ball size

            redCoorX = round(coorRect(3)/reduceFactor);
            redCoorY = round(coorRect(4)/reduceFactor);

            [x, y] = meshgrid(1:redCoorX,1:redCoorY);

            sizesU = unique(sizeV);

            spikeSum = zeros(size(Mr,1),length(respU),sizeX(4),'single');
            %spikeSumBase = mean(Mbase,3);

            %%%%%% Add spikes accordying to number of frames.
            msPerFarme= round(stimDur/sizeX(4));
            Delay = 250;
            [Mr] = BuildBurstMatrix(goodU,round(p.t),round((directimesSorted)),round((stimDur)));

            [trials neurons bins] = size(Mr);

            for u = 1:length(respU)

                Mu = squeeze(Mr(:,respU(u),:));

                j= 1;

                for f = 1:sizeX(4)

                    spikeSum(:,u,f) = sum(Mu(:,1*j:min(f*msPerFarme,length(Mu))),2);

                    j = f*msPerFarme;
                end
            end
            spikeSum = spikeSum/msPerFarme;
        end

        spikeSumArt = zeros(size(spikeSum));

        spikeSumArt(16:30,1,1:20) =1;

        %%%% Calculate grid positions %%%%%
        % Preallocate cell array to store indices of coordinates in each grid square

        ballRadius = sizesU(ceil(sizeN/2))/2; %Divided by for so that at least 1/4 of the ball needs to be inside the square  
        squarePresence = cell(gridSize, gridSize);
        [numTrials numPositions] = size(ChangePosX);

        % Loop through directions
        
        % Loop through trials
        for trial = 1:numTrials %(numTrials/direcN)*3+1:numTrials 
            for pos = 1:numPositions
                x = ChangePosX(trial, pos)-(screenSideO(3)-screenSideO(4))/2; % x-coordinate
                y = ChangePosY(trial, pos); % y-coordinate

                % Check which grid square(s) the ball overlaps
                for Yg = 1:gridSize
                    j=1;
                    for Xg = sort(1:gridSize,'descend')
                        % Get square boundaries
                        yMin = (Yg - 1) * squareSize; %- ballRadius;% - ballRadius;
                        yMax = Yg * squareSize; %- ballRadius;% + ballRadius;
                        xMin = (Xg - 1) * squareSize; %- ballRadius; %Half the ball has entered
                        xMax = (Xg) * squareSize; %- ballRadius;%the ball starts to leave.
                       

                        % Check if ball is within extended boundaries
                        if x >= xMin && x <= xMax && y >= yMin && y <= yMax
                            % Store the trial and time point for this square
                            squarePresence{Yg,j} = [squarePresence{Yg,j}; trial, pos];
                        end
                        j=j+1;
                    end
                end
            end
        end

        % Take mean of spike rate of each position

        spikeRatePos = zeros(length(respU),direcN,gridSize,gridSize);
        spikeSumArt = zeros(size(spikeSum));
        spikeSumArt(trialDivision*9*3+1:trialDivision*9*3+15,1,end-20:end) =1;%

        for u = 1:length(respU)

            for row = 1:gridSize
                for col = 1:gridSize
                    for d =4%1:direcN

                        spikeSumU = squeeze(spikeSum(C(:,2)==uDir(d),u,:));
                        spikeSumU = squeeze(spikeSumArt(C(:,2)==uDir(d),1,:));
                        rowIndices = squarePresence{row, col}(:,1);
                        rowIndices = rowIndices(ismember(rowIndices,find(C(:,2)==uDir(d))));
                        colIndices = squarePresence{row, col}(:,2);
                        colIndices = colIndices(ismember(rowIndices,find(C(:,2)==uDir(d))));
                        %Convert to linear indices
                        linearIndices = sub2ind(size(spikeSumU), ...
                            rowIndices-(trialDivision*offsetN)*(d-1), ...
                            colIndices);
                       % linearIndices = sub2ind(size(spikeSumU),rowIndices,colIndices);


                        spikeRatePos(u,d,row,col) = mean(spikeSumU(linearIndices));
                    end

                end

            end
        end
        %spikeRatePosP = permute(spikeRatePos,[1,2,4,3]).*1000;
        fig = tiledlayout(direcN/2,direcN/2);
        for d = 1:direcN
        nexttile;imagesc(squeeze(spikeRatePos(2,d,:,:)));colorbar;caxis([0 max(spikeRatePos(2,:,:,:),[],'all')]);
        title(string(uDir(d)))
 
        end
        figure;imagesc(squeeze(spikeSumArt(:,1,:)));%colormap(flipud(gray(64)));
        yline(trialDivision:trialDivision:size(spikeSumArt,1));
        yline(trialDivision*9:trialDivision*9:size(spikeSumArt,1)-1,'r','LineWidth',5)

    
        %%%%

        boot_means = load(sprintf('MovBall-Base-Boot-1000-%s',NP.recordingName)).boot_means;
        RFuSize =  load(sprintf('RFuSizeC-%s',NP.recordingName),'RFuSize').RFuSize;

        uSize

        ZScoreU = zeros(length(unique(seqMatrix)),size(goodU,2));


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

end
close all
end

