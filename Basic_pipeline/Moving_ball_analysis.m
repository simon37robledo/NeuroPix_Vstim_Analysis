cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%Optionall
summPlot = 0;
plotexamplesMB =0;
newTIC = 0;
ex=1;
ZscoresDo=1;
ReceptiveField = 0;

%%
% Iterate through experiments (insertions and animals) in excel file
for ex =1:size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
    %ex = 25;
    path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    cd(path)
    NP = NPAPRecording(path);

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
            else
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

    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,1,"MovBall",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,1,"MovBall",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A

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
    goodUdepth = NP.chLayoutPositions(2,goodU(1,:));

    %6. Get depths of units correcting for angle and micromanipulator depth:
    verticalDepth = sin(deg2rad(data.Angle(ex)))*(data.Depth(ex) - goodUdepth);

    %7. Load raster matrices
    bin = 1;
    preBase = 300;%round(3*interStimStats/4);

    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((stimDur+preBase*2)/bin)); %response matrix
    [nT,nN,nB] = size(Mr);

    %figure;imagesc(squeeze(Mr(:,15,:)));colormap(flipud(gray(64)));xline(preBase/bin);xline((preBase+stimDur)/bin);yline(trialDivision:trialDivision:nT-1);

    [MbTotal] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-interStimStats)/bin),round((interStimStats)/bin));% baseline matrix (whole baseline)
    baseline = squeeze(mean(squeeze(mean(MbTotal)),2));
    MrMean = squeeze(mean(Mr));

    v_replicated = repmat(baseline, 1, size(MrMean, 2));
    Norm = MrMean ./ v_replicated; %normalization by division
    Norm(Norm>2) = 2;

    %8. Plot summary response
    summPlot =0;
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
            if ~exist(path+"\Figs",'dir')
                mkdir Figs
            end
            cd(NP.recordingDir + "\Figs")
            print(fig, sprintf('%s-MovBall-summary.png',NP.recordingName),'-dpng');
        end
    end

    stims = stimOn';

    %%%%%%%%%%%%%% Select baseline


    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

    Mb = mean(Mb,3); %mean across time bins

    [trials, neurons] = size(Mb);

    % MbC = zeros(trials, neurons);

    %     percentTrials = round(trials*0.2);
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

    %spkRateBM = mean(MbC); %Calculate baseline spike rate of each neuron.

    %2. Calculate denominator for Z-score calculation
    spkRateBM = mean(Mb); %total mean.
    epsilon = 0.001;
    %denom = mad(MbC,0)+epsilon; %Calculate baseline variation of each neuron.

    denom = std(Mb)+eps;
    %%%%%%%%%%%%%%%%%%%% Select responsive units and calculate their tunning profile
    %3. Convolute matrix to spread response
    % Convolute in the 3rd dimension (trials)
if ZscoresDo ==1
    [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');

    [nT,nN,nB] = size(MrC);

    %4. Create window to scan rasters and get the maximum response
    duration =250; %ms
    window_size = [1, duration/bin];

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

    %Real data:
    for u =1:nN
        % Slide the window over the matrix
        %unit matrix
        max_mean_value(u) = -Inf; %General max? not needed
        NeuronRespProfile = zeros(nT/trialDivision,9); %4 columns plus: ofsett, dir, size, speed, frec.

        k =1;
        for i = 1:trialDivision:nT %Iterate across trials
            uM = mean(squeeze(MrC(i:i+trialDivision-1,u,:)));%Create matrix per unique trial conditions

            max_position_Trial = zeros(nT/trialDivision,2); %Create 2 matrices, for mean inside max window, and for position of window, for each trial category
            max_mean_value_Trial = zeros(1,nT/trialDivision);
            max_mean_value_Trial(k) = -Inf;

            for j = 1:size(uM, 2) - window_size(2) + 1 %Iterate across bins
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

            %4%. Z-score
            % NeuronRespProfile(i,4) = (max_mean_value_Trial(i) - (spkRateBM(u)+max_mean_value_Trial(i))/2)/denom(u); %Zscore
            NeuronRespProfile(k,4) = (max_mean_value_Trial(k) - spkRateBM(u))/denom(u); %Zscore
            %Assign visual stats to last columns of NeuronRespProfile. Select
            %according  to trial (d) the appropiate parameter > directions'offsets' sizes' speeds' freq'
            NeuronRespProfile(k,5) = C(i,2);
            NeuronRespProfile(k,6) = C(i,3);
            NeuronRespProfile(k,7) = C(i,4);
            NeuronRespProfile(k,8) = C(i,5);
            NeuronRespProfile(k,9) = C(i,6);

            k = k+1;

        end

        %figure;imagesc(uM);xline(max_position_Trial(i,2));xline(max_position_Trial(i,2)+window_size(2))

        NeuronVals(u,:,:) = NeuronRespProfile;
    end
    save(sprintf('NeuronRespCat-%s',NP.recordingName),"NeuronVals")

    %     [MrespO] = BuildBurstMatrix(goodU,round(p.t/bin),round(directimesSorted/bin),round((stimDur+duration)/bin));
    %     Z_scores = (MrespO-spkRateBM)./denom;
    %
    %     max(Z_scores,[],'all')

    %Create tunning curve, based on direction tunning curve --> clustering
    tuningCurve = zeros(nN,direcN);

    uDir = unique(directions);

    for u = 1:nN
        NeuronD = squeeze(NeuronVals(u,:,[4 5]));
        for d = 1:direcN
            tuningCurve(u,d) = NeuronD(NeuronD(:,2)==uDir(d),1)';
        end
    end

    save(sprintf('tuningC-%s',NP.recordingName),"tuningCurve")

    %Select significantly responsive units with shuffling plus Z-score:
    rands = 1000;

    if ~isfile(sprintf('randZscores-It-%d.mat',rands))
        RandZscoreU = zeros(rands, nN,2);
        for s = 1:rands
            %%%Shuffle trials
            M_shuffTr=Mr(randperm(size(Mr,1)),:,:); %check if neuron is responsive to specific stim (shuffle across trials)
            M_shufMTrC = ConvBurstMatrix(M_shuffTr,fspecial('gaussian',[1 5],3),'same');

            %%%Shuffle times
            ti2 = zeros(1,length(p.t));
            for u = 1:nN
                ti=p.t(p.ic(3,u):p.ic(4,u));
                diffi = diff(ti);
                %ti = ti(1) + cumsum(randperm(round(diff(ti))));
                ti2(p.ic(3,u):p.ic(4,u)) = [ti(1) ti(1)+cumsum((diffi(randperm(numel(diffi)))))];
            end

            M_shuffTi=BuildBurstMatrix(goodU,round(ti2/bin),round((stims-preBase)/bin),round((stimDur+preBase*2)/bin)); %check if neuron is visually response (shuffle across time)
            M_shufMTiC = ConvBurstMatrix(M_shuffTi,fspecial('gaussian',[1 5],3),'same');
            MB_shufTi = BuildBurstMatrix(goodU,round(ti2/bin),round((stims-preBase)/bin),round(preBase/bin));
            MB_shufTi = mean(MB_shufTi,3);
            spkRateBMS = mean(MB_shufTi); %total mean.
            %denom = mad(MbC,0)+epsilon; %Calculate baseline variation of each neuron.

            denomSti = std(MB_shufTi)+eps;

            NeuronRespProfileShTr = zeros(nT,nN);
            NeuronRespProfileShTi = zeros(nT,nN);
            for u =1:nN
                % Slide the window over the matrix (this selects the random
                % max)
                %unit matrix

                k =1;
                for i = 1:trialDivision:nT %Iterate across trials
                    uMTr = squeeze(M_shufMTrC(i+trialDivision-1,u,:))';%Create matrix per unique trial conditions
                    uMTi = squeeze(M_shufMTiC(i+trialDivision-1,u,:))';%Create matrix per unique trial conditions

                    %Create 2 matrices, for mean inside max window, and for position of window, for each trial
                    max_mean_value_Trialtr = zeros(1,nT);
                    max_mean_value_Trialtr(i) = -Inf;

                    %Create 2 matrices, for mean inside max window, and for position of window, for each trial
                    max_mean_value_Trialti = zeros(1,nT);
                    max_mean_value_Trialti(i) = -Inf;

                    for j = 1:size(uM, 2) - window_size(2) + 1 %Iterate across bins
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

                    %. Z-score
                    NeuronRespProfileShTr(i,u) = (max_mean_value_Trialtr(i) - spkRateBM(u))/denom(u);%Mb_shufTrM(u))/denomTr(u); %Zscore
                    NeuronRespProfileShTi(i,u) = (max_mean_value_Trialti(i) - spkRateBMS(u))/denomSti(u);%Mb_shufTiM(u))/denomTi(u); %Zscore

                    k=k+1;
                end

            end

            tDivList = 1:trialDivision:nT-1;

            meanTshuffledTr = zeros(u,numel(tDivList));
            meanTshuffledTi = zeros(u,numel(tDivList));

            %
            for td = 1:numel(tDivList) %take means across trial divisions to look at random times
                meanTshuffledTi(:,td) = mean(NeuronRespProfileShTi(tDivList(td):tDivList(td)+trialDivision-1,:));
                meanTshuffledTr(:,td) = mean(NeuronRespProfileShTr(tDivList(td):tDivList(td)+trialDivision-1,:));
            end

            %tri = randi(numel(tDivList)); %Take group of random trials to look at trials

            %RandZscore matrix: 3 dimension, 1st (number of iterations) 2nd
            %(number of neurons), 3rd (1. trials shuffled, 2. Times
            %shuffled).
            RandZscoreU(s,:,1) =  max(meanTshuffledTr,[],2);
            RandZscoreU(s,:,2) =  max(meanTshuffledTi,[],2);

            % RandZscoreU(s,:,2) =  mean(NeuronRespProfileShTi(tDivList(tri):tDivList(tri)+trialDivision-1,:));%max(meanTshuffledTi,[],2);
        end


        save(sprintf('randZscores-It-%d',rands),"RandZscoreU");
    else %if eand Z-score file exists then load

        RandZscoreU = load(sprintf('randZscores-It-%d.mat',rands));
        RandZscoreU = RandZscoreU.RandZscoreU;

    end

    Zthreshold = 3;
    maxZScores = squeeze(max(tuningCurve,[],2));
    %%Check if at least for 1 direction neuron is visually responsive:

    maxZS_vr = max(meanZScores,[],2);

    pvalTr = (sum(squeeze(RandZscoreU(:,:,1)) >= maxZScores') +1)/(rands+1);
    pvalTi = (sum(squeeze(RandZscoreU(:,:,2)) >= maxZScores') +1)/(rands+1);

    figure;histogram(squeeze(RandZscoreU(:,15,1)))

    save(sprintf('pvalTrials-%d',NP.recordingName),'pvalTr')
    save(sprintf('pvalTime-%d',NP.recordingName),'pvalTi')


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

    %Tuning index
    OI = sqrt(sum(meanZScoresSig.*sin(2*uDir),2).^2 + sum(meanZScoresSig.*cos(2*uDir),2).^2)./sum(meanZScoresSig,2); % OI calculated as in Dragoi and Swindale.

    a = sum()


%     nColors = 64; % Number of colors in the colormap
%     red = [linspace(0, 1, nColors)]'; % Red component increases linearly
%     green = [linspace(0, 0, nColors)]'; % Green component is always 0
%     blue = [linspace(1, 0, nColors)]'; % Blue component decreases linearly
%     customColormap = [red, green, blue]; % Combine into an nColors-by-3 matrix
%     normalizedValues = round(rescale([floor(min(verticalDepthSig')/100)*100; verticalDepthSig'; ceil(max(verticalDepthSig')/100)*100], 1, nColors)); %

    %Normalize and scale values to 1-64
    colorsForValues = customColormap(round(normalizedValues(2:end-1)), :); % Get corresponding color for each value

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

    %     % Assuming 'data' is your dataset
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


%     exampleU = NeuronVals{u};
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
end
    % Rasters plots per Neuron
    for plotOp = plotexamplesMB
        if plotexamplesMB == 1

            eNeuron =1:length(goodU);
            eNeuron = 1;

            orderS = [2 3 4 5;4 2 3 5;5 2 3 4;3 2 4 5];
            orderNames = {'dir_off_sizes_speeds';'sizes_dir_off_speeds';'speeds_dir_off_sizes';'off_dir_sizes_speeds'};

            A = [stimOn directions' offsets' sizes' speeds'];


            for k = 1

                [C sIndex2]= sortrows(A,orderS(k,:));

                %Sort directions:

                %directimesSorted = C(:,1)';

                [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((stimDur+preBase*2)/bin));

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
                    title(sprintf('U.%d-R.%.3f-B.%.3f-S.%.3f',u,max_mean_value(u),spkRateBM(u),Zscore(u)));

                    xticks([0.5 (preBase/bin):10:nB])
                    xticklabels([-preBase 0:10*bin:nB*bin])
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
                    rectangle('Position', [max_position(u,2)/mergeTrials, max_position(u,1)/mergeTrials, window_size(2)/mergeTrials, window_size(1)/mergeTrials],...
                        'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
                    hold off
                    prettify_plot

                    cd(NP.recordingDir)
                    if ~exist(path+"\Figs",'dir')
                        mkdir Figs
                    end
                    cd(NP.recordingDir + "\Figs")

                    fig.Position = [353    42   734   954];
                    print(fig, sprintf('%s-MovBall-%s-U%d-W%d-%dW.png',NP.recordingName,orderNames{k},u,window_size(1),window_size(2)),'-dpng');
                    %close

                end

            end
        end


    end
    %%%%%%%%%%%%%%%%%%%%%%%% Position heatmap

    if ReceptiveField ==1
        %%%

        %Get X and Y positions

        Xpos = cell2mat(ball.VSMetaData.allPropVal(22));

        Ypos = cell2mat(ball.VSMetaData.allPropVal(23));

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
            while tr <= size(ChangePosX,2)%length(MtrialCat) | 

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








%% Find shape of receptive field and whether it is signifficanlty different than the rest of the image. 
fun = @(block_struct) ...
   mean(block_struct.data,'all');
for u = 1:size(goodU,2)
    for s = 1:sizeN

        M=squeeze(normRFu(s,:,:,u));
        M2 = blockproc(M,[5 5],fun);
        NaNs = isnan(M);
        M(isnan(M))=-1;
        M2 = blockproc(M,[5 5],fun);
        %M=M+1;

        figure;imagesc(M);colorbar
        M2n = M;%normalize(M2,'center','mean');

        maxBlock = max(M2n(~NaNs),[],'all');
        minBlock = min(M2n(~NaNs),[],'all');

        respDelta = sum(maxBlock-M2n((~NaNs)),'all',"omitnan")/numel(M2n((~NaNs)));
        RTI = 1-((respDelta)^2/(maxBlock-minBlock)^2);

    end

end


%%

% Load an image
img = squeeze(normRFu(:,:,17));

% Convert to grayscale if it's not
if size(img,3) == 3
    img = rgb2gray(img);
end

% Define block size
blockSize = [50, 50]; % example size

% Define the function to find maximum
fun = @(block_struct) max(block_struct.data(:));

% Apply the function to each block
maxBlock = blockproc(img, blockSize, fun);

% Display the result
imshow(maxBlock, []);

figure(normRFu(normRFu~=-1))
 td= squeeze(normRFu);
 td = td(td~=-1);
 figure;imagesc(td);

figure;imagesc(squeeze(normRFu(:,:,:,17)));colorbar
figure;imagesc(squeeze(RFu(:,:,:,17)));colorbar
figure;imagesc(squeeze(matrixResp(:,:,20)))
figure;imagesc(squeeze(matrixForNorm));colorbar
figure;imagesc(matrix_nXnY);

unique(matrix_nXnY)
unique(matrixResp)
unique(RFu(:,:,20))

colorbar


%% Calculate the spike sum per frame

spikeSum = zeros(nT,length(goodU),sizeX(4));
%spikeSumBase = mean(Mbase,3);

% Make raster matrix that contains only stimulus. 
msPerFarme= round(stimDur/sizeX(4));
[Mr] = BuildBurstMatrix(goodU,round(p.t),round((directimesSorted)),round((stimDur)));
[Mbase] = BuildBurstMatrix(goodU,round(p.t),round((directimesSorted-preBase)),round((preBase)));

for u = 1:length(goodU)

    Mu = squeeze(Mr(:,u,:));

    j= 1;  

    for f = 1:sizeX(4)  
        
        spikeSum(:,u,f) = sum(Mu(:,1*j:min(f*msPerFarme,length(Mu))),2);

        j = f*msPerFarme;
    end
end
spikeSum = single(spikeSum);
spikeRates=spikeSum./(msPerFarme/1000); %Divide by seconds per frame
    

%% Prepare for convolution
%%get videos of every trial
% Define the grid size
xAxis = 1:cell2mat(ball.VSMetaData.allPropVal(63));%63 for pv103,64 for pv67 %65 PV139 and on   %:max(X,[],'all')-min(X,[],'all')-max(sizeV)/2;

yAxis = 1:cell2mat(ball.VSMetaData.allPropVal(63));%63 PV103,64 for pv67 %65 and on   %max(Y,[],'all')-min(Y,[],'all')-max(sizeV)/2;

%respU = [9,11,12,13,15,16,17,18,19,20,21,22,23,24,25,31]; %PV139_2

%respU =[18,21,22,24,26,30,31,32,34,35,38,39,41,42,44,47,50,51,53,57,58,61,63,65,66]; %PV139_1

%respU = [66, 65, 61, 42, 24,  21]; %PV59_1

%respU = [24]; %PV32_1
%respU = [2,8,16]; %PV32_2
%respU = [45,9]; %PV32_3
% respU = [18,22,24,26,42,46,47,50,53,56,61,64,65,68,73,79,81,84,89,90,94,97,98,100,101,102,103,104,107,108,112,...
%      113,119,120,121,124,127,130,131,133,137,138,141,142,144,145,147,149,151,154,157]; %PV67_1
%respU = [115,102,77,75]; %PV67_4
respU = [1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 17, 20, 22, 23, 24, 28, 32, 34, 36]; %PV103_1
%respU = [2,23,25,28,20,22,25,28,29,32,34,36,37,38,39,40]; %PV103_2
%respU = [3,7,10,14,15,17,18,19,20,26,33,46,50]; %PV103_3
%respU = [5,8,13];%PV103_4
%respU = [35,89];%PV103_5
%respU = [11,13,14,21,29,31,32,28,45,80,81,82]; %PV103_6
%respU = [6,18,25,29,30,32,36,37,38,44,57,78,82,84,89]; %PV103_7
%respU = sort(respU);

cd(NP.recordingDir + "\Figs")

fixed_delay = 1;


RFuDir = cell(1,direcN);
RFuSize = cell(1,sizeN);
RFuSpeed = cell(1,speedN);

A = [stimOn directions' offsets' sizes' speeds'];

C = sortrows(A,[5 2 3 4]);

uDir = unique(directions);
uSize = unique(sizes);
uSpeed = unique(speeds);

for i = 1:direcN
    C(C(:,2) == uDir(i),2) = i;
end
for i = 1:sizeN
    C(C(:,4) == uSize(i),4) = i;
end
for i = 1:speedN
    C(C(:,5) == uSpeed(i),5) = i;
end

centerXreal = cell2mat(ball.VSMetaData.allPropVal(61)); %PV139 63 %PV64 62 %PV103 61 
%% Convolution
%respU =[144,145,147,149,151,154,157];


%respU = 18;
for u = 1:length(respU)
    
    RFu = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');

    for k = 1:direcN
        RFuDir{k} = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');
    end
    for k = 1:sizeN
        RFuSize{k} = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');
    end
    for k = 1:speedN
       RFuSpeed{k} = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');
    end

    for i = 1:trialDivision:trials

        videoTrials = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');
        for j = 1:sizeX(4) %%Calculate video of unique trials

            % Create a matrix of zeros representing the grid
            xyScreen = zeros(length(xAxis)+1,length(yAxis)+1,'single');

            %Define the center and radius of the circle

%             centerX = ChangePosX(i,j)-min(X,[],'all'); % X position of the center
%             centerY = ChangePosY(i,j)-min(Y,[],'all'); % Y position of the center
            centerX = ChangePosX(i,j)-(centerXreal-length(yAxis)/2); % X position of the center
            centerY = ChangePosY(i,j); % Y position of the center

            radius = sizeV(i)/2; % Radius of the circle

            % Create a meshgrid of coordinates
            [x, y] = meshgrid(1:length(xAxis)+1, 1:length(yAxis)+1);

            % Calculate the distance of each point from the center
            distances = sqrt((x - centerX).^2 + (y - centerY).^2);

            % Set the values inside the circle to 1 (or any other value you prefer)
            xyScreen(distances <= radius) = 1;

            videoTrials(:,:,j) = xyScreen;
            %figure;imagesc(xyScreen);
        end

        spikeMean = mean(spikeSum(i:i+trialDivision-1,respU(u),:))-mean(Mbase(i:i+trialDivision-1,respU(u),:),'all'); %Get the mean across same trials and substract the baseline

        tic
        Co = convn(videoTrials,spikeMean,'same');
        toc
        %toc

%         spikeMean = mean(spikeSum(i:i+trialDivision-1,respU(u),:));
% 
%         tic
%         Co = convn(videoTrials,spikeMean,'same');
%         toc

        RFuDir{C(i,2)} =  RFuDir{C(i,2)}+Co;
        RFuSize{C(i,4)} = RFuSize{C(i,4)}+Co;
        RFuSpeed{C(i,5)} = RFuSpeed{C(i,5)}+Co;

        RFu = RFu+Co;
    end

    save(sprintf('RFu-%d',respU(u)),'RFu')
    save(sprintf('RFuDir-%d',respU(u)),'RFuDir','-v7.3')
    save(sprintf('RFuSize-%d',respU(u)),'RFuSize','-v7.3')
    save(sprintf('RFuSpeed-%d',respU(u)),'RFuSpeed','-v7.3')
end
%%




%%
%Delay = 204 ms = 12 bins of 17 ms

delay = 12;
if rem(sizeX(4),2)==0
    center = sizeX(4)/2+1;
else
    center = (sizeX(4)-1)/2+1;
end

%%Search for maximal intensity
% figure
% plot([1:154],squeeze(spikeMean)')

%% Depth distribution of responses

basicPathA = {basic_pathPV67,basic_pathPV103,basic_pathPV27,basic_pathPV139,basic_pathPV59,basic_pathPV32};

expA = {expPV67,expPV103,expPVPV27,expPV139,expPV59,expPV32};

InsAnglesPV67 = [88 88 88 88];
InsDepthsPV67 = [3956, 3907, 4000, 3934]; 
RespUnitPV67 ={[18,22,24,26,42,46,47,50,53,56,61,64,65,68,73,79,81,84,89,90,94,97,98,100,101,102,103,104,107,108,112,...
    113,119,120,121,124,127,130,131,133,137,138,141,142,144,145,147,149,151,154,157],[],[],[115,102,77,75]};

InsAnglesPV103= [70.5 70.5 70.5 70.5 70.5 70.5 70.5 ];
InsDepthsPV103= [4104.14, 3964.18, 4066.5, 4123.87, 4175.86, 4225.55, 4027.42]; 
RespUnitPV103={[1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 17, 20, 22, 23, 24, 28, 32, 34, 36],[2,23,25,28,20,22,25,28,29,32,34,36,37,38,39,40],...
    [3,7,10,14,15,17,18,19,20,26,33,46,50],[5,8,13],[35,89],[11,13,14,21,29,31,32,28,45,80,81,82],[6,18,25,29,30,32,36,37,38,44,57,78,82,84,89]};

InsAnglesPV27= [72.5 72.5 72.5  72.5 72.5 72.5 72.5];
InsDepthsPV27= [3913.9 3904.34 3525.7 3900.1 3914.34 3739.3 3906.8];
RespUnitPV27={};

InsAnglesPV139 = [89 89];
InsDepthsPV139 = [3910 3907.23];
RespUnitPV139 = {[18,21,22,24,26,30,31,32,34,35,38,39,41,42,44,47,50,51,53,57,58,61,63,65,66],...
    [9,11,12,13,15,16,17,18,19,20,21,22,23,24,25,31]};

InsAnglesPV59 = [81 81 81];
InsDepthsPV59 = [2845.89 3358.91 3050.36];
RespUnitPV59 = {[66, 65, 61, 42, 24,  21],[],[]};

InsAnglesPV32 = [69 72 72];
InsDepthsPV32 = [2400 2600 2200];
RespUnitPV32 = {[24],[2,8,16],[45]};

animalD = {InsDepthsPV67,InsDepthsPV103,InsDepthsPV27,InsDepthsPV139,InsDepthsPV59,InsDepthsPV32};

animalA = {InsAnglesPV67,InsAnglesPV103,InsAnglesPV27,InsAnglesPV139,InsAnglesPV59,InsAnglesPV32};

RespUnits = {RespUnitPV67,[RespUnitPV103],[],RespUnitPV139,RespUnitPV59,RespUnitPV32};

verticalDepth = cell(2,length(animalA));

YDist = cell(2,length(animalA));
%% Depth plot
figure()
j = 1;
aj=2;
for a = [1,2,4:length(animalA)] %X doesn't change


    %point1 = [animals{a}(:,1),animals{a}(:,2),animals{a}(:,3)];

    verticalDepth{1,a} = sin(deg2rad(animalA{a})).*(animalD{a}); %depth of unit along vertical axis

    YDist{1,a} = cos(deg2rad(animalA{a})).*(animalD{a}); %X distance of unit from insertion

    %point2 = [animals{a}(:,1),animals{a}(:,2)+YDist',animals{a}(:,3)-verticalDepth'];

    %     for i = 1:length(animals{a})
    %         plot3([point1(i,1),point2(i,1)],[point1(i,2),point2(i,2)],[point1(i,3),point2(i,3)],'LineWidth', 1, 'Color',colorA{a})
    %         hold on
    %     end

    for in = 1:length(animalA{a})
        path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
            in+string(filesep)+"catgt_"+string(expA{a})+"_"+in+"_g0");

        NP = NPAPRecording(path);

        cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

        %Good units

        GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

        GoodU_orDepth = cluster_info.depth(cluster_info.group=="good");

        %GoodU_orDepth = GoodU_orDepth(RespUnits{a}{in});

        verticalDepth{2,a}{in} = sin(deg2rad(animalA{a}(in)))*(animalD{a}(in) - GoodU_orDepth);

        ResponsiveDepth =  verticalDepth{2,a}{in}(RespUnits{a}{in});
        Nindex = true(1, length(verticalDepth{2,a}{in}));
        Nindex(RespUnits{a}{in}) = false;
        NonRespDepth =  verticalDepth{2,a}{in}(Nindex);

        hold on;
        plot(ones(length(ResponsiveDepth),1)+j+0.5-1+aj-2, ResponsiveDepth,'o','MarkerFaceColor','b')
        hold on; plot(ones(length(NonRespDepth),1)+j-1+aj-2, NonRespDepth,'o','MarkerFaceColor','r')

        j=j+3;
    end
aj = aj+5;
end

ylabel('Unit depth (um)')
xlabel('insertions and animals (blue= visually responive units to Moving ball stim)')


%% find unique pairs and their counts. Create Normalization matrix
%coords = [ChangePosX(:)-min(X,[],'all')+1,ChangePosY(:)-min(Y,[],'all')+1];
%centerXreal = cell2mat(ball.VSMetaData.allPropVal(61)); %63 in PV139 and on
%coords = [ChangePosX(:)+1-(centerXreal-length(yAxis)/2),ChangePosY(:)+1];%translate coordinates to 1080x1080

%task: ask Mark about exact parameters of screen, what is the appropiate parameter? 

difX = max(ChangePosX,[],'all') - min(ChangePosX,[],'all'); %max(X,[],'all')-min(X,[],'all')
difY = max(ChangePosY,[],'all') - min(ChangePosY,[],'all');
% 
% nX = max(coords(:,1));
% nY = max(coords(:,2));
%[Cn, ~, ic] = unique(coords, 'rows', 'stable');
%counts = accumarray(ic, 1);
find(arrayfun(@(x) strcmp(x.allPropName , 'rect'), ball.VSMetaData));
coorRect = cell2mat(ball.VSMetaData.allPropVal(57));
find(contains(ball.VSMetaData.allPropName,'rect'));%63 in PV139 and on

%[x, y] = meshgrid(1:length(xAxis)+1, 1:length(yAxis)+1);
[x, y] = meshgrid(1:coorRect(3),1:coorRect(4));
%[x, y] = meshgrid(1:difX,1:difY);
sizesU = unique(sizeV);

matrixForNorm = zeros(length(xAxis)+1, length(yAxis)+1);
%matrixForNorm  = zeros(1920+1,1920+1);
%matrixForNorm  = zeros(difX,difY);

matrixForNorm = zeros();

for s = 1:sizeN
    % Step 3: Fill in the matrix with the counts of repeated coordinates
    for i = 1:trialDivision*sizeN:length(ChangePosX)
        %matrix_nXnY = zeros(length(xAxis)+1, length(yAxis)+1);
        matrix_nXnY = zeros(coorRect(4),coorRect(3));
        %matrix_nXnY = zeros(difX,difY);
        for tr = 1:sizeX(4)

            %             centerX = Cn(i,1);
            %             centerY = Cn(i,2);

            centerX = ChangePosX(i,tr);%-(centerXreal-length(yAxis)/2);
            %centerX = ChangePosX(i,tr)*length(xAxis)/coorRect(3);
            centerY = ChangePosY(i,tr);
            %centerY = (ChangePosY(i,tr)+(centerXreal-length(yAxis)/2))*length(yAxis)/coorRect(3);

            % Create a meshgrid of coordinate
            radius = sizesU(s)/2;

            % Calculate the distance of each point from the center
            distances = sqrt((x - centerX).^2 + (y - centerY).^2);

            % Set the values inside the circle to 1 (or any other value you prefer)
            matrix_nXnY(distances <= radius) =1;
           
        end
        matrixForNorm = matrixForNorm+matrix_nXnY;
    end
      
end

figure;imagesc(matrixForNorm)
colorbar
%needs to be transformed into a 1080x180
%PV103 done / PV64 / PV139 
expA = {expPV67,expPV103,expPVPV27,expPV139,expPV59,expPV32};
cd(NP.recordingDir)
save(sprintf('matrixForNorm-%s',expA{4}),'matrixForNorm')


%%%Message for mark:
% Regarding the moving ball graphs, when I try to reconstruct a summary image of the stimulus, I get something weird. Either not all offsets are shown, or all offsets are there but some don't intersect, as if the extreme offsets where off the screen (I checked by running a stim with the same parameters, and indeed all offsets intersect). This is a problem of identifying on which coordinate reference the X and Y of the moving ball are recorded. I am using  the rect parameter which is a [0 0 1920 1080]. This results in offsets not intersecting. 
%   


%% Plot RF with a stationary delay
%to do:
%1. Standarize  the screen size to 1080 (for animal 6, bigger screen size)
%2. Detect maximum screen value per unit, and normalize, then sum across
%units. 
%3. Find better way to detect delay in convolution (max of video?). 
%4. Save result matrices. 


cd(NP.recordingDir+"\Figs")
%respU = 25;
% respU
% %respU = 157;
sideVF = 1080;%cell2mat(ball.VSMetaData.allPropVal(63)); %65 in pv139 and on

% SumAllDir = cell(2,length(animalA)); 
% SumPerDir = cell(2,length(animalA)); 
sizeScreenAn = [64,63,64,65,65,65]; %changes in the parameters order. 
%%
for a =[2]
    inSumAll = cell(1,length(animalA{a}));
    
    inSumDir = cell(2,length(animalA{a}));
    

    for in = 1:length(animalA{a})
        respU = RespUnits{a}{in};

        %generate path
        path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
            in+string(filesep)+"catgt_"+string(expA{a})+"_"+in+"_g0");
        NP = NPAPRecording(path);
        
        %Load stim stats
        patternIndex = strfind(string(NP.recordingDir), "\catgt");
        endIndex = patternIndex(1)-1;
        stimDir = string(NP.recordingDir);
        stimDir = extractBetween(stimDir,1,endIndex);
        file = dir (stimDir);
        filenames = {file.name};
        ballFiles = filenames(contains(filenames,"linearlyMovingBall"));
        ball= load(stimDir+"\"+string(ballFiles));
        directions = cell2mat(ball.VSMetaData.allPropVal(17));
        uDirections = unique(directions);
        sizeScreen = cell2mat(ball.VSMetaData.allPropVal(sizeScreenAn(a)));
        
        cd(path+"\Figs")

        Allu = zeros(sizeScreen+1,sizeScreen+1,length(respU));
        AlluDir = zeros(sizeScreen+1,sizeScreen+1,length(respU),length(uDirections));        

        for u =1:length(respU)
            
            %respU(u) =157;
            %sideVF = cell2mat(ball.VSMetaData.allPropVal(63)); %65 in pv139 and on

            M = load(sprintf('RFu-%d.mat',respU(u)));
            M = M.RFu;
            Mdir = load(sprintf('RFuDir-%d.mat',respU(u)));
            Mdir= Mdir.RFuDir;

            [sM1 sM2 sM3] = size(M);

            if sM1 >  sideVF+1 || sM2 >  sideVF+1

                leftPix = floor((sM1-sideVF)/2);

                rightPix = sM1-floor((sM1-sideVF)/2)-1;

                M = M(leftPix:rightPix,leftPix:rightPix,:);

            end

%             fig = figure; % Open a new figure window
%             imagesc(squeeze(M(:,:,center+delay)))%.
%             %axis equal; % Ensure the aspect ratio is equal to make the circle look r ound
%             %axis off; % Optionally, turn off the axis for a cleaner look
%             title(sprintf('U.%d-ReceptiveField-Del-%d ms',respU(u),delay*msPerFarme));
%             set(gca,'YDir','normal')
%             colorbar
%             %prettify_plot
%             cd(path+"\Figs")
%             %print(fig, sprintf('%s-MovBall-U%d-RFconv.png',NP.recordingName,respU(u)),'-dpng');
%             close

            Allu(:,:,u) = M(:,:,center+delay);

          %  t = tiledlayout(length(Mdir)/2,length(Mdir)/2);

            for direc = 1:length(Mdir)

                %nexttile; % Open a new figure window
                Md = Mdir{direc};
               % imagesc(squeeze(Md(:,:,center+delay)))%.
                %axis equal; % Ensure the aspect ratio is equal to make the circle look r ound
                %axis off; % Optionally, turn off the axis for a cleaner look
                %title(sprintf('U.%d-ReceptiveField-Del-%d ms',respU(u),delay*msPerFarme));
                %set(gca,'YDir','normal')
                %colorbar
                %caxis([0, 0.6]);
                %prettify_plot

                AlluDir(:,:,u,direc) =  Md(:,:,center+delay);
                %dirsSum{2,direc} = uDirections(direc);
                
            end     
%             cd(path+"\Figs")
%             print(fig, sprintf('%s-MovBall-U%d-RFconv.png',NP.recordingName,respU(u)),'-dpng');
%             close

        end
        inSumAll{in} =  Allu;
        inSumDir{1,in} = AlluDir;
        inSumDir{2,in} = uDirections;
    end
    SumAllDir{1,a} = inSumAll;
    SumAllDir{2,a} = expA{a};
    SumPerDir{1,a}= inSumDir;
    SumPerDir{2,a}= expA{a};
end


%% Get summary images. 


for a = [1,2,4]

    f1= figure;
    subplot(1,3,1)

    for in=1:length(animalA{a})

        uA = SumAllDir{1,a}{in};
        meanA = mean(uA,3);

        %imagesc(meanA-


    end

end
%%



M = load('RFu-26.mat');
M = M.RFu;
[aaa, aaaain] = max(M,[],3);

[aa,aain] = max(aaa,[],'all');

aaaain(aain); 
fig = figure; % Open a new figure window
imagesc(squeeze(M(:,:,center+delay)))%./matrixForNorm*spkRateBM(respU));
%colormap('gray'); % Use a grayscale colormap





M = load('RFu26.mat');
M = M.RFu;

% %% Test conv center
% spikeMean(1,1,:) = zeros(1,1,154); 
% spikeMean(1,1,78) = 1; %--> center = 78 154/2+1
% C = convn(videoTrials,spikeMean,'same');



% %%Video File Test
% %M = rand(240, 320, 100);
% videoFile = 'u29.avi';
% writerObj = VideoWriter(videoFile);
% writerObj.FrameRate = 60; % Set the frame rate to 30 frames per second
% open(writerObj);
% 
% for k = 1:size(M, 3) % Loop through each frame
%     frame = M(:, :, k); % Extract the k-th frame
%     % If your data isn't uint8, you might need to normalize and convert it
%     %frame = im2uint8(mat2gray(frame)); % Normalize and convert to uint8
%     writeVideo(writerObj, frame); % Write the frame to the video
% end
% 
%    close(writerObj);
% 
%    implay(videoFile);

%%
    
% %%Iterate per neuron
% 
% for u = 1:length(goodU)
% 
%     Mu = squeeze(Mr(:,u,:));
% 
%     j= 1;
% 
%     spikeSum = zeros(trials,sizeX(4));
% 
%     for f = 1:sizeX(4)  
%         
%         spikeSum(:,f) = sum(Mu(:,1*j:min(f*msPerFarme,length(Mu))),2);
% 
%         j = f*msPerFarme;
% 
%     end
% 
%     iter = max(spikeSum,[],'all'); %iteration number
% 
%     for i = 1:iter
% 
%         [row,col] = find(spikeSum>=i);
% 
%         Ysum = ChangePosY(row,col);
%         Xsum = ChangePosX(row,col);
% 
%         X_flat = Xsum(:); % Flatten A into a column vector
%         Y_flat = Ysum(:); % Flatten B into a column vector
% 
%         % Combine into pairs
%         pairs = [X_flat, Y_flat];
% 
%         % Find unique pairs and their counts
%         [unique_pairs, ~, ic] = unique(pairs, 'rows', 'stable');
%         pair_counts = accumarray(ic, 1);
% 
%         %Transoform pair intro screen coordinates
%         x = round(unique_pairs(:,1))-round(min(X,[],'all'))+1;
%         y = round(unique_pairs(:,2))-round(min(Y,[],'all'))+1;
% 
%         for c = 1:length(pair_counts)
%             xyScreen(y(c), x(c)) = pair_counts(c);
%         end
% 
%         nGrid = 5;
%         
%         redGrid = zeros(nGrid,nGrid);
% 
%         for b = 1:nGrid
% 
%             for v = 1:nGrid
% 
%                 div = round(length(xyScreen)/nGrid);
% 
%                 redGrid(b,v) = mean(xyScreen(b*div-div+1:b*div,v*div-div+1:v*div),'all');
% 
%             end
% 
%         end
%         figure
%         
%         imagesc(redGrid);
% 
%         max(y)
% 
% 
%     end
% 
% 
%     for x = 1:length(xAxis)
%         for y = 1:lrngth(yAxis)
% 
%             xpos=find
% 
%         end
%     end
% 
% 
% 
% end
% 


