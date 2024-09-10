cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%Optionall
summPlot = 0;
plotexamplesMB =0;
newTIC = 0;
ZscoresDo=1; redoResp=1;
Shuffling =1;
repeatShuff =1;
ReceptiveFieldFixedDelay = 0;
tuning =1;
depthPlot =0;
ReceptiveFieldConvolutions =0;
x=1;
examplesSDG =[28 29 30 31 40 41 42 43]; %[1 2 3 4 5 6 7 8 9 10 11 12 13 14 29 30 31 32 40 41 42 43]:
summPlot =0;

%%%In shuffling make sure that response cat is selected equally between SDG
%%%and MB
%%
% Iterate through experiments (insertions and animals) in excel file
for ex = examplesSDG%examplesSDG%[7 8 28]%1:size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
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
    ballFiles = filenames(contains(filenames,"linearlyMovingBall"));

    
        if isempty(ballFiles)
            %disp()
            w= sprintf('No moving ball files where found in %s. Skipping into next experiment.',NP.recordingName);
            warning(w)
            continue
        end

    directions = [];
    offsets = [];
    sizes = [];
    speeds = [];
    orientations = [];

    j =1;


    if size(ballFiles) ~= [0 0]

        for i = ballFiles %Extract stim properties
            ball= load(stimDir+"\"+string(i));

            if ~isempty(find(strcmp(ball.VSMetaData.allPropName,'orientations'))) %Check if orientations are present (grid inside moving object).
                orientations = [orientations cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'orientations'))))];
            else %No orientations property exist
                orientations = [orientations zeros(1,cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials')))))];
            end
            if isempty(orientations) %orientations property exists but is empty
                orientations = [orientations zeros(1,cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials')))))];
            end

            directions = [directions cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'directions'))))];

            offsets = [offsets cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'offsets'))))];

            sizes = [sizes cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballSizes'))))];

            speeds = [speeds cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'speeds'))))];

            interStimStats = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    if ~isempty(find(strcmp(ball.VSMetaData.allPropName,'orientationFreq')))
        Freq = [ cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'orientationFreq'))))];
    end

    
    uDir = unique(directions);
    offsetN = length(unique(offsets));
    direcN = length(unique(directions));
    speedN = length(unique(speeds));
    sizeN = length(unique(sizes));
    orientN = length(unique(orientations));
    nT = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials'))));
    trialDivision = nT/(offsetN*direcN*speedN*sizeN*orientN); %Number of trials per unique conditions
    categ = nT/trialDivision;


    %3. Load Triggers (diode)
    Ordered_stims= strsplit(data.VS_ordered{ex},',');
    containsMB = cellfun(@(x) contains(x,'MB'),Ordered_stims);
    ttlInd = find(containsMB);

    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,1,"MB",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,1,"MB",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A

    %Check diode
    % for i = 1:length(stimOn)
    %
    %     Lon = length(onSync(onSync>=stimOn(i) & onSync<=stimOff(i)));
    %
    %     Loff = length(offSync(offSync>=stimOn(i) & offSync<=stimOff(i)));
    %
    %     if Lon ~= 77 || Loff ~= 77
    %
    %         2+2
    %
    %     end
    %
    %     a = NP.getAnalogData(1,stimOn(i)-500,4000);
    %
    %     figure;plot(squeeze(a));
    %     xline((onSync(onSync>=stimOn(i) & onSync<=stimOff(i))-stimOn(i)-500)*NP.samplingFrequencyNI/1000)
    %
    % end

    stimInter= mean(stimOn(2:end)-stimOff(1:end-1)); % When dealing with different speeds, save different stim durs, and create M for each speed
    stimDur = mean(-stimOn+stimOff); %When dealing with different speeds, save different stim durs, and create M for each speed
    A = [stimOn directions' offsets' sizes' speeds' orientations'];
    [C indexS] = sortrows(A,[2 3 4 5 6]);

    %4. Sort directions:
    directimesSorted = C(:,1)';

    B = [stimOff directions' offsets' sizes' speeds' orientations'];
    [Coff indexSo] = sortrows(B,[2 3 4 5 6]);
    directimesSortedOff = Coff(:,1)';

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

    p = NP.convertPhySorting2tIc(NP.recordingDir);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');
    goodSPKS = p.nSpks(label == 'good');
    goodAmp = p.neuronAmp(label == 'good');
    %goodUdepth = NP.chLayoutPositions(2,goodU(1,:));

    %6. Get depths of units correcting for angle and micromanipulator depth:
    %verticalDepth = sin(deg2rad(data.Angle(ex)))*(data.Depth(ex) - goodUdepth);

    %7. Load raster matrices
    bin = 1;
    preBase = round(stimInter/2);%round(3*interStimStats/4);

    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-stimInter/2)/bin),round((stimDur+stimInter)/bin)); %response matrix
    [nT,nN,nB] = size(Mr);

    %figure;imagesc(squeeze(Mr(:,15,:)));colormap(flipud(gray(64)));xline(preBase/bin);xline((preBase+stimDur)/bin);yline(trialDivision:trialDivision:nT-1);

    [MbTotal] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-interStimStats)/bin),round((interStimStats)/bin));% baseline matrix (whole baseline)
    baseline = squeeze(mean(squeeze(mean(MbTotal)),2));
    MrMean = squeeze(mean(Mr));

    v_replicated = repmat(baseline, 1, size(MrMean, 2));
    Norm = MrMean ./ v_replicated; %normalization by division
    Norm(Norm>2) = 2;

    %8. Plot summary response

    for plotOp = summPlot
        if summPlot

            fig = figure;
            imagesc(Norm);
            xline(preBase/bin,'k', LineWidth=1.5)
            xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5)
            hcb = colorbar();
            title(hcb,'SpkR/Baseline');
            xticks([0.5 (preBase/bin):10:nB])
            xticklabels([-preBase 0:10*bin:nB*bin])
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
            cd(NP.recordingDir)
            cd(NP.recordingDir + "\Figs")
            print(fig, sprintf('%s-MovBall-summary.png',NP.recordingName),'-dpng');
        end
    end

    stims = stimOn';

    %%%%%%%%%%%%%% Select baseline
    NP.recordingDuration_ms/1000/60

    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

    Mb = mean(Mb,3); %mean across time bins

%     percentTrials = round(nT*0.2);
% 
%     [S, Sindex] = sort(Mb);
% 
%     % Mb = S(1:end-percentTrials,:);
% 
%     [trials, neurons] = size(Mb);
% 
%     %randomize again.
%     %Mb = Mb(randi([1, size(Mb,1)], 1, size(Mb,1)),:);
% 
%     %Take mean across window of % of trials
% 
%     MbC = zeros(round(trials/percentTrials), neurons);
% 
%     for i = 1:percentTrials:trials
%         meanb = mean(Mb(i:min(i+percentTrials-1, end),:),1);
%         MbC(j,:) = meanb;
%         j = j+1;
%     end
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
    for ZscoreData =1
if ZscoresDo ==1
    Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((stimDur)/bin)); %response matrix
    [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');

    [nT,nN,nB] = size(MrC);

    %4. Create window to scan rasters and get the maximum response
    duration =400; %ms
    window_size = [1, round(duration/bin)];

    %5. Initialize the maximum mean value and its position

    max_position = zeros(nN,2);
    max_mean_value = zeros(1,nN);

    %     for u =1:nN %Initial window
    %         % Slide the window over the matrix
    %         %unit matrix
    %         max_mean_value(u) = -Inf;
    %         for d = 1:length(unique(offsets))*trialDivision:nT
    %             uM = squeeze(MrC(d:d+length(unique(offsets))*trialDivision-1,u,:));
    %             for i = 1:size(uM, 1) - window_size(1) + 1
    %                 for j = 1:size(uM, 2) - window_size(2) + 1
    %                     % Extract the sub-matrix
    %                     sub_matrix = uM(i:min(i+window_size(1)-1, end), j:min(j+window_size(2)-1,end));
    %
    %                     % Compute the mean value
    %                     mean_value = mean(sub_matrix(:));
    %
    %                     % Update the maximum mean value and its position
    %                     if mean_value > max_mean_value(u)
    %                         max_mean_value(u) = mean_value;
    %                         max_position(u,:) = [i+d-1, j];
    %
    %                     end
    %
    %                 end
    %             end
    %         end
    %     end

    NeuronVals = zeros(nN,nT/trialDivision,9); %Each neuron, which has a matrix where the first column is maxVal of bins, 2nd and 3rd position of window in matrix...
    % 4th Z-score.
    % responsive compared to the baseline.

    %[1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 17, 20, 22, 23, 24, 28, 32, 34, 36]

    if ~isfile(sprintf(sprintf('NeuronRespCat-%s.mat',NP.recordingName))) || redoResp
    
    %Baseline = size window

    [Mbd] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-duration)/bin),round(duration/bin)); %Baseline matrix plus

    %Merge trials:

    mergeTrials = trialDivision;
% 
%     Bd = reshape(Mbd, [mergeTrials, size(Mbd, 1)/mergeTrials, size(Mbd, 2), size(Mbd,3)]);
% 
%     Mbd2 = mean(squeeze(mean(Bd, 1)),3);


    B = reshape(MrC, [mergeTrials, size(Mr, 1)/mergeTrials, size(Mr, 2), size(Mr,3)]);

    % Take the mean across the first dimension (rows)
    Mr2 = squeeze(mean(B, 1));

    [nT,nN,nB] = size(Mr2);
    cd(NP.recordingDir)

    %Real data:
    tic
    for u =1:nN
        % Slide the window over the matrix
        %unit matrix
        max_mean_value(u) = -Inf; %General max? not needed
        NeuronRespProfile = zeros(nT,9); %4 columns plus: ofsett, dir, size, speed, frec.

        k =1;
        max_position_Trial = zeros(nT,2); %Create 2 matrices, for mean inside max window, and for position of window, for each trial category
        max_mean_value_Trial = zeros(1,nT);
        for i = 1:nT %Iterate across trials
            uM = squeeze(Mr2(i,u,:))';%Create matrix per unique trial conditions
            %uMb = Mbd2(i,u);
           
            max_mean_value_Trial(k) = -Inf;

            for j = 1:2:size(uM, 2) - window_size(2) + 1 %Iterate across bins
                % Extract the sub-matrix
                sub_matrix = uM(j:min(j+window_size(2)-1,end)); %Create matrix of size window per bin

                % Compute the mean value
                mean_value = mean(sub_matrix(:)); %Compute mean of each window

                % Update the maximum mean value and its position (a
                % window is selected across each trial division
                if mean_value >  max_mean_value_Trial(k)
                    max_mean_value_Trial(k) = mean_value;
                    max_position_Trial(k,:) = [i j];
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
            NeuronRespProfile(k,4) = max_mean_value_Trial(k)-spkRateBM(u);
            %Assign visual stats to last columns of NeuronRespProfile. Select
            %according  to trial (d) the appropiate parameter > directions'offsets' sizes' speeds' freq'
            NeuronRespProfile(k,5) = C(i*mergeTrials,2);
            NeuronRespProfile(k,6) = C(i*mergeTrials,3);
            NeuronRespProfile(k,7) = C(i*mergeTrials,4);
            NeuronRespProfile(k,8) = C(i*mergeTrials,5);
            NeuronRespProfile(k,9) = C(i*mergeTrials,6);

            k = k+1;

        end

        %figure;imagesc(uM);xline(max_position_Trial(i,2));xline(max_position_Trial(i,2)+window_size(2))

        NeuronVals(u,:,:) = NeuronRespProfile;
    end
        save(sprintf('NeuronRespCat-%s',NP.recordingName),"NeuronVals")
    else
        NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;
    end
toc
    %Create tunning curve, based on direction tunning curve 
    tuningCurve = zeros(nN,direcN);

    udir = unique(directions);

    for u = 1:nN
        NeuronD = squeeze(NeuronVals(u,:,[1 5]));
        for d = 1:direcN
            tuningCurve(u,d) = max(NeuronD(NeuronD(:,2)==udir(d),1))'; %Selecting top direction. 
        end
    end

    save(sprintf('tuningC-%s',NP.recordingName),"tuningCurve")
%     unnn = 1:67;
%     unnn(any(tuningCurve > 4, 2))


%Keep real max window for shuffling.
[mp mi] = max(NeuronVals(:,:,1),[],2);

psT = zeros(size(mi, 1), 1);
psB = zeros(size(mi, 1), 1);

% Extract the values from A(:,:,2:3) according to the max indices
for i = 1:size(mi, 1)
    psT(i) = NeuronVals(i, mi(i), 2);
    psB(i) = NeuronVals(i, mi(i), 3);
end
    
for Shuffle =1
    if Shuffling ==1
        rands = 1000;
        [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-stimInter/2)/bin),...
            round((stimDur+stimInter)/bin)); %response matrix
        sMr = single(Mr);
        %sMb = single(Mbd);
%         [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round(preBase/bin));
%         sMb = single(Mb);
        if ~isfile(sprintf('randValues-It-%d.mat',rands))||repeatShuff==1
            %Select significantly responsive units with shuffling plus Z-score:
           tic
           numSegments = 4;
            parfor par = 1:numSegments
                RandValU = zeros(rands/numSegments, nN,2,'single');
                tic

                for s = 1:rands/numSegments
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
                    
                    %%%Shuffle baseline
%                     B = reshape(sMb(randperm(size(sMb,1)),:,randperm(size(sMb,3))), [mergeTrials, size(sMb, 1)/mergeTrials, size(sMb, 2), size(sMb,3)]);
%                     MB_shuff =  squeeze(mean(B, 1));
%                     MB_shuff = squeeze(mean(MB_shuff,3));
                    
                    NeuronRespProfileShTr = zeros(nN,'single');
                    NeuronRespProfileShTi = zeros(nN,'single');
                    NeuronRespProfileShDiff = zeros(nN,'single');


                    for u =1:nN
                        %                     Slide the window over the matrix (this selects the random
                        %                     max)
                        %                     unit matrix

                        k =1;
                        for i = 1:nT %Iterate across trials
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

                RandValUp{par} = RandValU;

            end
            toc
            RandValU = cat(1, RandValUp{:});
            
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
end
end
    end


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

for tuningIndex = 1
if tuning ==1

    cd(NP.recordingDir)
    %pvalTr=  load(sprintf('pvalTrials-%s',NP.recordingName)).pvalTr;
    %pvalTi= load(sprintf('pvalTime-%s',NP.recordingName)).pvalTi;
    tuningCurve = load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve;
    NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;

    %Tuning index
    % Calculate tuning direction and strength for all neurons then plot it across depth and in 3D
    uDir = unique(rad2deg(directions));
    uDirRad = deg2rad(uDir);
    OI = sqrt(sum(tuningCurve.*sin(2*uDirRad),2).^2 + sum(tuningCurve.*cos(2*uDirRad),2).^2)./sum(tuningCurve,2); % OI calculated as in Dragoi and Swindale.

    %Preferred orientation
    oneSideCurve = zeros(nN,direcN/2);
   

    %     for u = 1:nN
    %         NeuronD = squeeze(NeuronVals(u,:,[1 5]));
    %         NeuronD(:,2)=rad2deg(NeuronD(:,2));
    %         for d = 1:length(uDir)/2
    %             oneSideCurve(u,d) = mean(NeuronD(NeuronD(:,2)==uDir(d),1))';
    %             oneSideCurve(u,d) = oneSideCurve(u,d)+ mean(NeuronD(NeuronD(:,2)==uDir(d)+180,1))';
    %         end
    %     end

    a = sum(tuningCurve.*cos(2*uDirRad),2);
    b = sum(tuningCurve.*sin(2*uDirRad),2);

    Theta = zeros(1,nN);
    for u = 1:nN
        if a(u)>0
            Theta(u) = rad2deg(0.5*atan(b(u)/a(u)));
        else
            Theta(u) = rad2deg(deg2rad(180)+0.5*atan(b(u)/a(u)));
        end
    end

    L = sqrt(a.^2+b.^2)./sum(tuningCurve,2); %Tuning Strength

    [preferDir pI] = max(tuningCurve,[],2);

    nonPreferDir = zeros(nN,1);
    for u = 1:nN
        npI = find(uDir == mod(uDir(pI(u)) + 180, 360));
        nonPreferDir(u) = tuningCurve(u,npI);
    end
    DSI = 1- nonPreferDir./preferDir;
    DSI(isnan(DSI))=0;
    L(isnan(L))=0;

    if any(DSI>1)

        2+2

    end
    

    save(sprintf('Angle-prefer-%s',NP.recordingName),'Theta')
    save(sprintf('Orientation-Tuning-Index-%s',NP.recordingName),'L')
    save(sprintf('Direction-Selectivity-Index-%s',NP.recordingName),'DSI')



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

%Clustering Analysis:
%     BinZScore = meanZScoresSig >Zthreshold;
%
%     %normS = normalize(meanZScoresSig);
%     %normS = diff(log10(meanZScoresSig+1+abs(min(meanZScoresSig,[],'all'))));
%     normS = log(meanZScoresSig+1+abs(min(meanZScoresSig,[],'all')));
%     %figure;imagesc(squeeze(Mr(:,1,:)));colormap(flipud(gray(64)))
%
%     %Clustering tuning curves:
%     %
%     %     k = 4; % Example number of clusters
%     %     [idx, Cl] = kmeans(BinZScore,
%     %     k,'MaxIter',10000,'Display','final','Replicates',10,'OnlinePhase','on');
%
%
%    [idx] = dbscan(normS,1,2); %seems to work better than k means, no need to know number of clusters
%
%
%
%     %idx = clusterdata(meanZScoresSig,2);
%     %figure;plot(rad2deg(uDir), normS, 'LineWidth', 2);
%
%     k = unique(idx);
%
%     % Visualize the clusters (assuming 2D data for simplicity)
%     fig = figure;
%     gscatter(normS(:,1), normS(:,2), idx);
%     title('Neuron Clustering with K-means');
%     xlabel('Feature 1');
%     ylabel('Feature 2');
%     legend('Cluster 1', 'Cluster 2', 'Cluster 3');
%     cd('\\132.66.45.127\data\Large_scale_mapping_NP\Immobilized_exp\Orientation_tuning_figs');
%     print(fig, sprintf('%s-MovBall-clusters_groups.png',NP.recordingName),'-dpng');
%     close all
%
%     %Plot depth as color across different clusters
%
%     fig = figure;
%     for i = 1:length(k)
%         subplot(1,length(k),i)
%         pl= plot(rad2deg(uDir), meanZScoresSig(idx ==k(i),:), 'LineWidth', 2);
%         xlabel('Angle (degrees)');
%         xticks(rad2deg(uDir));
%         xticklabels(arrayfun(@num2str, rad2deg(uDir), 'UniformOutput', false));
%         ylabel('Mean Z-Score Response');
%         title(sprintf('Tuning Curve cluster %d',i));
%         grid on;
%         ylim([floor(min(meanZScoresSig,[],'all')) ceil(max(meanZScoresSig,[],'all'))])
%
%         clusterColors = colorsForValues(find(idx ==k(i)),:);
%
%         for j = 1:length(pl)
%             set(pl(j),'color',clusterColors(j,:)); % Change colour of plot.
%         end
%
%     end
%
%
%     colormap(flipud(customColormap))
%     cb = colorbar();
%     cb.TickLabels = flip(linspace(floor(min(verticalDepthSig)/100)*100,ceil(max(verticalDepthSig)/100)*100,11));
%     cb.Label.String = 'Depth';
%     cb.Label.FontSize = 12;
%     cd('\\132.66.45.127\data\Large_scale_mapping_NP\Immobilized_exp\Orientation_tuning_figs');
%     print(fig, sprintf('%s-MovBall-clustersR.png',NP.recordingName),'-dpng');
%close all

% Enhance the plot with error bars showing the standard error of the mean
%     sem = stdResponses / sqrt(numTrials); % Standard Error of the Mean
%     hold on;
%     errorbar(angles, meanZScores, sem, 'LineStyle', 'none', 'Color', 'k');
%     hold off;
%
%     sse = [];
%     for k = 1:10
%         [idx, C, sumd] = kmeans(meanZScores, k);
%         sse(end+1) = sum(sumd);
%     end
%     figure;
%     plot(1:10, sse, '-o');
%     xlabel('Number of clusters k');
%     ylabel('Sum of squared distances');
%     close all

% Assuming 'data' is your dataset
%     k = 3; % Example number of clusters
%     [idx, C] = kmeans(meanZScores, k);
%
%     % Visualize the clusters (assuming 2D data for simplicity)
%     figure;
%     gscatter(meanZScores(:,1), meanZScores(:,2), idx);
%     title('Neuron Clustering with K-means');
%     xlabel('Feature 1');
%     ylabel('Feature 2');
%     legend('Cluster 1', 'Cluster 2', 'Cluster 3');

% exampleU = NeuronVals{u};
%
%     spkRateR = max_mean_value;
%
%     mSpk = mean(spkRateR);
%
%     Zscore = (spkRateR - (spkRateBM + spkRateR)/2)./denom;
%
%     RespU = GoodU_or(Zscore>3);
%     RespUInd =  find(Zscore>3);
%
%     sDirections = sort(directions);
%
%     %    tunning = sDirections(max_position(:,1));
%
%     MinusRB = spkRateR - spkRateBM;

%K-means clustering --> Look at sidekick rply
%1% preprocessing the data
%- Pass each neuron into a row? Use 3D matrix directly? Or use raster plot as image? Ask Mark.
%- Normalize the data (Z-score?)
%-

% Rasters plots per Neuron
for plotOp = plotexamplesMB
    if plotexamplesMB == 1

        eNeuron =34;%1:length(goodU); %8
        %eNeuron = 1;

        orderS = [2 3 4 5;4 2 3 5;5 2 3 4;3 2 4 5];
        orderNames = {'dir_off_sizes_speeds';'sizes_dir_off_speeds';'speeds_dir_off_sizes';'off_dir_sizes_speeds'};

        A = [stimOn directions' offsets' sizes' speeds'];

        cd(NP.recordingDir)
        NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;

        posX = squeeze(NeuronVals(:,:,3));
        posY = squeeze(NeuronVals(:,:,2));

        %2%. WindowTrial
        %             NeuronRespProfile(k,2) = max_position_Trial(k,1);
        %
        %             %3%. WindowBin
        %             NeuronRespProfile(k,3) = max_position_Trial(k,2);

        uDir = unique(rad2deg(directions));
        bin0 = bin;
        bin =50;
        trialsPerAngle = trialDivision*offsetN*speedN*sizeN*orientN;
        for k = 1

            [C sIndex2]= sortrows(A,orderS(k,:));

            %Sort directions:

            directimesSorted = C(:,1)';

            [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((stimDur+preBase*2)/bin));

            Mr2 = [];

            for u = eNeuron

                j=1;

                mergeTrials = 1;

                if mergeTrials ~= 1 %Merge trials

                    for i = 1:mergeTrials:trials

                        meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end),u,:)),1);

                        Mr2(j,:) = meanb;

                        j = j+1;

                    end
                else
                    Mr2 = Mr(:,u,:);
                end

                [nT,nB] = size(Mr2);

                fig = figure;

                imagesc(squeeze(Mr2));colormap(flipud(gray(64)));
                %Plot stim start:
                xline(preBase/bin,'k', LineWidth=1.5)
                %Plot stim end:
                xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5)
                ylabel('Trials');xlabel('Time (ms)');
                title(sprintf('U.%d-Unit-phy-%d',u,GoodU_or(u)));

                %xticks([0.5 (preBase/bin):10:nB])
                %xticklabels([-preBase 0:10*bin:nB*bin])
                yticklabels([yticks]*mergeTrials)
                %Directions
                v = nT/direcN:nT/direcN:nT-1;
                yline(v+0.5,'r', LineWidth=3);
                %Offsets
                v = nT/(direcN*offsetN):nT/(direcN*offsetN):nT-1;
                yline(v+0.5,'b', LineWidth=2);
                %sizes
                v = nT/(direcN*offsetN*sizeN):nT/(direcN*offsetN*sizeN):nT-1;
                yline(v+0.5, LineWidth=0.5);

                %                             hcb = colorbar();
                %                             title(hcb,'Spikes/sec');
                %caxis([0 max(0.2,max(max_mean_value(u)))])
                hold on
                %Plot rectangle:
                for d = 1:size(NeuronVals,2)
                rectangle('Position', [posX(u,d)/(bin/bin0)+round(preBase/bin), (posY(u,d)*trialsPerAngle-trialsPerAngle)/mergeTrials, window_size(2)/(bin/bin0), trialDivision],...
                    'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
                end
                hold off
                prettify_plot
                
                yyaxis right
                ylim([1,nT])
                yticks([trialsPerAngle:trialsPerAngle:nT])
                ax = gca;
                ax.YAxis(2).Color = 'k';
                yticklabels(sort(uDir,'descend'))
                ylabel('Angles')
                lims =xlim;
                xt = xticks;

                cd(NP.recordingDir)
                if ~exist(path+"\Figs",'dir')
                    mkdir Figs
                end
                cd(NP.recordingDir + "\Figs")

                fig.Position = [353    42   734   954];
                print(fig, sprintf('%s-MovBall-%s-U%d-W%d-%dW.png',NP.recordingName,orderNames{k},u,window_size(1),window_size(2)),'-dpng');
                close

                %%%Polar plot
                tuningCurve = (load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve)*1000; %convert to spikes/sec
                theta = deg2rad(uDir); %linspace(0, 2*pi, size(tuningCurve,2)+1);  % 9 points for 8 bars (because it's circular)
                  % Remove the last value to avoid duplication of the first
                pf = figure;
                % Create the polar plot
                polarplot([theta, theta(1)], [tuningCurve(eNeuron,:), tuningCurve(eNeuron,1)], '-o')
                set(pf,"Color",'w')
                ax = gca;
                ax.ThetaTick = uDir;
                title(sprintf('PolarPlot-U.%d-Unit-phy-%d',u,GoodU_or(u)));
                print(fig, sprintf('%s-MovBall-polarPlot-U%d.png',NP.recordingName,u),'-dpng');
                close


            end

        end
    end


end
%%%%%%%%%%%%%%%%%%%%%%%% Position heatmap

for receptiveField=1 %%ReceptiveFieldMethod1

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

            [Mb] = BuildBurstMatrix(goodU,round(p.t),round((stims-preBase)),round(preBase)); %Baseline matrix plus

            Mb = mean(Mb,3); %mean across time bins

            spkRateBM = mean(Mb);

            A = [stimOn directions' offsets' speeds' orientations', sizes'];
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
                [MtrialCat] = BuildBurstMatrix(goodU,round(p.t),round((directimesSorted(i:i+trialDivision*sizeN-1))+delayResp),...
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

    verticalDepth = sin(deg2rad(data.Angle(ex))).*(data.Depth(ex)); %depth of unit along vertical axis

    YDist = cos(deg2rad(data.Angle(ex))).*(data.Depth(ex)); %X distance of unit from insertion

    verticalDepthU= verticalDepth- goodUdepth;

    NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName));
    tuningCurve = load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve;

    Theta = load(sprintf('Angle-prefer-%s',NP.recordingName)).Theta;
    TI= load(sprintf('Tuning-Index-%s',NP.recordingName)).L;

    %figure;polarplot(deg2rad(Theta()),L,'.','MarkerSize',20);set(gcf, 'Color', 'w')

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

%%%%% Convolution
for convNeuron = 1
    if ReceptiveFieldConvolutions ==1

        if ~isfile(sprintf('RFuC-%s.mat',NP.recordingName))%exist
            %%A. CONVOLUTION SETUP
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

            [x, y] = meshgrid(1:coorRect(3)/10,1:coorRect(4)/10);

            sizesU = unique(sizeV);

            %%%%%%%%%%Load responsive neurons
            pvalTi= load(sprintf('pvalTime-%s',NP.recordingName)).pvalTi;

            respUi = find(pvalTi <0.02);
            respU = [];

            for u = respUi

                percentT = sum(any(squeeze(Mr(:,u,:)>0),2))/size(Mr,1);

                if percentT >0.03 %Eliminate neurons that don't spike often.
                    respU = [respU u];
                end

            end


            spikeSum = zeros(size(Mr,1),length(respU),sizeX(4),'single');
            %spikeSumBase = mean(Mbase,3);

            %%%%%% Add spikes accordying to number of frames.
            msPerFarme= round(stimDur/sizeX(4));
            [Mr] = BuildBurstMatrix(goodU,round(p.t),round((directimesSorted)),round((stimDur)));
            [Mbase] = BuildBurstMatrix(goodU,round(p.t),round((directimesSorted-preBase)),round((preBase)));

            [trials neurons bins] = size(Mr);

            for u = 1:length(respU)

                Mu = squeeze(Mr(:,respU(u),:));

                j= 1;

                for f = 1:sizeX(4)

                    spikeSum(:,u,f) = sum(Mu(:,1*j:min(f*msPerFarme,length(Mu))),2);

                    j = f*msPerFarme;
                end
            end
            spikeSum = spikeSum/msPerFarme; %convert into spike rate.

            [Mb] = BuildBurstMatrix(goodU(:,respU),round(p.t),round((stims-preBase)),round(preBase)); %Baseline matrix plus

            Mb = mean(Mb,3); %mean across time bins

            spkRateBM = mean(Mb);

            % % %%artificial mod
            % delay= 10;
            % spikeSum = spikeSum(:,1,:);
            % spikeSum(1:15,1,31+delay:32+delay) = ones(15,1,2)+300;
            %
            % spikeSum(241+600:600+241+14,1,178-32+delay:178-31+delay) = ones(15,1,2)+300;

            %figure;imagesc(squeeze(spikeSum(:,1,:)))

            %respU =1;

            %%%%Initialize 5D matrices
            RFu = zeros(coorRect(4)/10,coorRect(3)/10,sizeX(4),length(respU),"single");

            RFuSpeed = zeros(speedN,coorRect(4)/10,coorRect(3)/10,sizeX(4),length(respU),"single");
            RFuDir = zeros(direcN,coorRect(4)/10,coorRect(3)/10,sizeX(4),length(respU),"single");
            RFuSize = zeros(sizeN,coorRect(4)/10,coorRect(3)/10,sizeX(4),length(respU),"single");


            Uspeed = unique(speeds);
            Udir = unique(directions);
            Usize = unique(sizes);
            NormVideo = zeros(coorRect(4)/10,coorRect(3)/10,sizeX(4),'single');

            %%A. CONVOLUTION. Runs the spike train across the stimulus videos, to
            %%extract noisy receptive field (because spike responses are noise).
            tic
            for i = 1:trialDivision:trials

                videoTrials = zeros(coorRect(4)/10,coorRect(3)/10,sizeX(4),'single');
                for j = 1:sizeX(4) %%Calculate video of unique trials


                    xyScreen = zeros(coorRect(4)/10,coorRect(3)/10,"single");
                    %matrixResp = zeros(coorRect(4),coorRect(3),nN,"single");

                    centerX = ChangePosX(i,j)/10;
                    centerY = ChangePosY(i,j)/10;
                    radius = sizeV(i)/2;

                    % Calculate the distance of each point from the center
                    distances = sqrt((x - centerX).^2 + (y - centerY).^2);

                    % Set the values inside the circle to 1 (or any other value you prefer)
                    xyScreen(distances <= radius/10) = 1;

                    videoTrials(:,:,j) = xyScreen;
                    %figure;imagesc(xyScreen);
                end

                spikeMean = mean(spikeSum(i:i+trialDivision-1,:,:)-spkRateBM); %Get the mean across same trials and substract the baseline
                Co = zeros(coorRect(4)/10,coorRect(3)/10,sizeX(4),length(respU),'single');


                for u = 1:length(respU)
                    Co(:,:,:,u) = convn(videoTrials,spikeMean(:,u,:),'same');
                end


                %     figure;imagesc(squeeze(Co(:,:,80,3)))
                RFuDir(Udir == C(i,2),:,:,:,:) = squeeze(RFuDir(Udir == C(i,2),:,:,:,:))+Co;
                RFuSize(Usize == C(i,6),:,:,:,:) = squeeze(RFuSize(Usize == C(i,6),:,:,:,:))+Co;
                RFuSpeed(Uspeed == C(i,4),:,:,:,:) = squeeze(RFuSpeed(Uspeed == C(i,4),:,:,:,:))+Co;

                NormVideo = NormVideo+videoTrials;

                RFu = RFu+Co;
            end
            toc
            save(sprintf('RFuC-%s',NP.recordingName),'RFu','-v7.3')
            save(sprintf('RFuDirC-%s',NP.recordingName),'RFuDir','-v7.3')
            save(sprintf('RFuSizeC-%s',NP.recordingName),'RFuSize','-v7.3')
            save(sprintf('RFuSpeedC-%s',NP.recordingName),'RFuSpeed','-v7.3')

        else

            RFu = load(sprintf('RFuC-%s',NP.recordingName),'RFu').RFu;
            RFuDir = load(sprintf('RFuDirC-%s',NP.recordingName),'RFuDir').RFuDir;
            RFuSize =  load(sprintf('RFuSizeC-%s',NP.recordingName),'RFuSize').RFuSize;
            RFuSpeed = load(sprintf('RFuSpeedC-%s',NP.recordingName),'RFuSpeed').RFuSpeed;

        end

        end
    end
end



% 
% set(gca, 'YDir', 'reverse');
% xticks([1:5:5*13])
% xticklabels(1:13)
% grid on
% set(gcf, 'Color', 'w')


% 
% implay(squeeze((RFuSizeSav(1,:,:,:,13))))
% 
% implay(squeeze(convolved_result_reshaped(:,:,:,13)))
%% %%%%%%% Find max receptive field and how localized it is (entropy): It runs
%%%%%%% a ball shaped window that moves across the image of the previously
%%%%%%% calculated convolution. 
RFuSizeSav = load(sprintf('RFuSizeC-%s',NP.recordingName)).RFuSize;
RFuDirSav = load(sprintf('RFuDirC-%s',NP.recordingName)).RFuDir; 

%Size analysis -> 1. hacer la grafica que pide Mark, diferentes Neuronas. 2. Calcular entropia en experimento y graficar contra la prfundidad. 3. Escribirle a sebastian caldas. 

convMeanSize = zeros(sizeN,coorRect(4)/10,coorRect(3)/10,length(respU));
maxValdelayS = zeros(sizeN,length(respU));
maxInddelayS = zeros(sizeN,length(respU));

for s=1:sizeN

    % Parameters for the circular mask
    radius = Usize(s)/10; % Define the radius of the circular window
    % Create a circular mask
    [X, Y] = meshgrid(-radius:radius, -radius:radius);
    circular_mask = (X.^2 + Y.^2) <= radius^2;

    % Normalize the circular mask so that the sum of the mask elements is 1
    circular_mask = double(circular_mask);
    circular_mask = circular_mask / sum(circular_mask(:));

    % Use convolution to find the sum of elements inside the circular window
    A = squeeze(RFuSizeSav(s,:,:,:,:));
    A_reshaped = reshape(A, size(A, 1), size(A, 2), []);

    %convolved_result = convn(A, reshape(circular_mask,[size(circular_mask,1),size(circular_mask, 2), 1, 1]), 'same');

    convolved_result = convn(A_reshaped, circular_mask, 'same');
    % Reshape the convolved result back to the original 4D structure
    convolved_result_reshaped = reshape(convolved_result, size(A, 1), size(A, 2), size(A, 3), size(A, 4));

    padding =150;
    A_padded = padarray(A, [0, 0, padding,0], 'post');

    % Note: The filter should be expanded to 3D by adding a singleton dimension
    circular_mask_3d = reshape(circular_mask, [size(circular_mask, 1), size(circular_mask, 2), 1, 1]);

    % Perform the convolution
    convolved_result_reshaped = convn(A_padded, circular_mask, 'same');


    % Initialize arrays to store results
    max_means = zeros(size(convolved_result, 3), size(convolved_result, 4));
    center_row = zeros(size(convolved_result, 3), size(convolved_result, 4));
    center_col = zeros(size(convolved_result, 3), size(convolved_result, 4));

    % Find the maximum values and their positions for each slice across the 3rd and 4th dimensions
    for j = 1:size(convolved_result_reshaped, 4)
        for i = 1:size(convolved_result_reshaped, 3)
            % Extract the 2D slice
            slice = convolved_result_reshaped(:, :, i, j);
            % Find the maximum value and its position
            [max_means(i, j), max_idx] = max(slice(:));
            [max_row, max_col] = ind2sub(size(slice), max_idx);
            % Adjust positions to get the center of the circular mask
            center_row(i, j) = max_row + radius;
            center_col(i, j) = max_col + radius;
        end
    end

    [max_across_3rd, idx_3rd] = max(max_means, [], 1);

    maxValdelayS(s,:) = max_across_3rd;

    maxInddelayS(s,:) = idx_3rd;

    for u = 1:length(respU)
        convMeanSize(s,:,:,u) =  squeeze(convolved_result_reshaped(:,:,idx_3rd(u),u));
    end
end
A_normalized = squeeze(convMeanSize(1,:,:,:));
for k = 1:length(respU)
    slice = A_normalized(:, :, k);
    min_val = min(slice(:));
    max_val = max(slice(:));
    A_normalized(:, :, k) = (slice - min_val) / (max_val - min_val);
end

[m i] = max(max_means(:,27))
figure;plot(max_means(:,27)');xline(90);xline(i,'Color','r');xlabel('frames');ylabel('max mean value')
implay(A_normalized)

%% Dir analysis


maxValdelayD = zeros(direcN,length(respU));
maxInddelayD = zeros(direcN,length(respU));

for d=1:direcN

    % Parameters for the circular mask
    radius = Usize(s)/10; % Define the radius of the circular window
    % Create a circular mask
    [X, Y] = meshgrid(-radius:radius, -radius:radius);
    circular_mask = (X.^2 + Y.^2) <= radius^2;

    % Normalize the circular mask so that the sum of the mask elements is 1
    circular_mask = double(circular_mask);
    circular_mask = circular_mask / sum(circular_mask(:));

    % Use convolution to find the sum of elements inside the circular window
    A = squeeze(RFuSizeSav(s,:,:,:,:));
    A_reshaped = reshape(A, size(A, 1), size(A, 2), []);

    %convolved_result = convn(A, reshape(circular_mask,[size(circular_mask,1),size(circular_mask, 2), 1, 1]), 'same');

    convolved_result = convn(A_reshaped, circular_mask, 'same');
    % Reshape the convolved result back to the original 4D structure
    convolved_result_reshaped = reshape(convolved_result, size(A, 1), size(A, 2), size(A, 3), size(A, 4));

    % Initialize arrays to store results
    max_means = zeros(size(A, 3), size(A, 4));
    center_row = zeros(size(A, 3), size(A, 4));
    center_col = zeros(size(A, 3), size(A, 4));

    % Find the maximum values and their positions for each slice across the 3rd and 4th dimensions
    for j = 1:size(A, 4)
        for i = 1:size(A, 3)
            % Extract the 2D slice
            slice = convolved_result_reshaped(:, :, i, j);
            % Find the maximum value and its position
            [max_means(i, j), max_idx] = max(slice(:));
            [max_row, max_col] = ind2sub(size(slice), max_idx);
            % Adjust positions to get the center of the circular mask
            center_row(i, j) = max_row + radius;
            center_col(i, j) = max_col + radius;
        end
    end

    [max_across_3rd, idx_3rd] = max(max_means, [], 1);

    maxValdelayS(s,:) = max_across_3rd;

    maxInddelayS(s,:) = idx_3rd;

end




%%


%Size
sizeSPK = zeros(sizeN,size(Mr,1)/sizeN,length(respU),sizeX(4));

dirsSPK =  zeros(direcN,size(Mr,1)/direcN,length(respU),sizeX(4));

for s=1:sizeN
    sizeSPK(s,:,:,:) = spikeSum(Usize(s)==C(:,6),:,:);
end

%Dir
for d = 1:direcN
 dirsSPK(d,:,:,:) = spikeSum(Udir(d)==C(:,2),:,:);
end

figure;imagesc(squeeze(sizeSPK1(:,27,:)));colormap(flipud(gray(64)));

nTs = size(spikeSum,1)/sizeN;
nT = size(spikeSum,1);
nTred = size(spikeSum,1)/trialDivision;

MeanbyTD = squeeze(mean(reshape(spikeSum,[trialDivision, size(spikeSum,1)/trialDivision,size(spikeSum,2),size(spikeSum,3)])));

MeanbyTDsize = squeeze(mean(reshape(sizeSPK,[trialDivision, size(sizeSPK,1), size(sizeSPK,2)/trialDivision,size(sizeSPK,3),size(sizeSPK,4)])));

figure;imagesc(squeeze(sizeSPK(1,:,27,:)));
colormap(flipud(gray(64)));
%Plot stim start:
%xline(preBase/bin,'k', LineWidth=1.5)
%Plot stim end:
%xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5)
ylabel('Trials');xlabel('Time (ms)');
%yticklabels([yticks]*mergeTrials)
%Directions
v = nTs/direcN:nTs/direcN:nTs-1;
yline(v+0.5,'r', LineWidth=3);
%Offsets
v = nTs/(direcN*offsetN):nTs/(direcN*offsetN):nTs-1;
yline(v+0.5,'b', LineWidth=2);
% %sizes
% v = nT/(direcN*offsetN*sizeN):nT/(direcN*offsetN*sizeN):nT-1;
% yline(v+0.5, LineWidth=0.5);



figure;imagesc(squeeze(spikeSum(:,27,:)));
colormap(flipud(gray(64)));
%Plot stim start:
%xline(preBase/bin,'k', LineWidth=1.5)
%Plot stim end:
%xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5)
ylabel('Trials');xlabel('Time (ms)');
yticklabels([yticks]*mergeTrials)
%Directions
v = nT/direcN:nT/direcN:nT-1;
yline(v+0.5,'r', LineWidth=3);
%Offsets
v = nT/(direcN*offsetN):nT/(direcN*offsetN):nT-1;
yline(v+0.5,'b', LineWidth=2);
%sizes
v = nT/(direcN*offsetN*sizeN):nT/(direcN*offsetN*sizeN):nT-1;
yline(v+0.5, LineWidth=0.5);


figure;imagesc(squeeze(MeanbyTD(:,27,:)));
colormap(flipud(gray(64)));
%Plot stim start:
%xline(preBase/bin,'k', LineWidth=1.5)
%Plot stim end:
%xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5)
ylabel('Trials');xlabel('Time (ms)');
yticklabels([yticks]*mergeTrials)
%Directions
v = nTred/direcN:nTred/direcN:nTred-1;
yline(v+0.5,'r', LineWidth=3);
%Offsets
v = nTred/(direcN*offsetN):nTred/(direcN*offsetN):nTred-1;
yline(v+0.5,'b', LineWidth=2);
%sizes
v = nTred/(direcN*offsetN*sizeN):nTred/(direcN*offsetN*sizeN):nTred-1;
yline(v+0.5, LineWidth=0.5);

%figure;imagesc(reshape(spikeSum,[trialDivision, size(spikeSum,1)/trialDivision,size(spikeSum,2),size(spikeSum,1)]))

[m in] = max(squeeze((RFuSize(1,:,:,:,27))),[],'all')
min(squeeze((RFuSize(1,:,:,:,27))),[],'all')
implay(squeeze((RFuSize(1,:,:,:,27))))

e = zeros(length(respU),sizeX(4));
for i =1:sizeX(4)
    for u = 1:length(respU)
        e(u,i) = entropy(double(squeeze(RFuSize(1,:,:,i,u))));
    end
end
    
size(RFuSize)

figure;plot(e(20,:))

[m in]= max(e(10,:))

min(squeeze(RFuSize(1,:,:,:,20)),[],'all')

