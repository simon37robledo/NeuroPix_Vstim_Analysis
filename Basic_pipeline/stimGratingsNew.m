%% stimGratings

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

newTIC = 0;
repeatShuff =0;
%%
%SA5_1,PV103_7,PV27_1

for ex =[7 8 28]%1:size(data,1)

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

    if ~exist(path+"\Figs",'dir')
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
    stimFiles = filenames(contains(filenames,"StaticDrifting"));

    
        if isempty(stimFiles)
            %disp()
            w= sprintf('No static- drifting gratings ball files where found in %s. Skipping into next experiment.',NP.recordingName);
            warning(w)
            continue
        end

    directions = [];
    tempFR = [];
    spatFR = [];

%     stim = load(stimDir+"\"+string(i));
%         static_time = cell2mat(grat.VSMetaData.allPropVal(11))*1000; %Static time
%         angles = [angles cell2mat(grat.VSMetaData.allPropVal(26))]; %Angles
% 
%         tf = [tf cell2mat(grat.VSMetaData.allPropVal(27))]; %time Freq
%         sp = [sp cell2mat(grat.VSMetaData.allPropVal(28))]; %spatial Freq
%         interStimStats = cell2mat(grat.VSMetaData.allPropVal(40))*1000;
    


    if size(stimFiles) ~= [0 0]

        for i = stimFiles %Extract stim properties
            stim= load(stimDir+"\"+string(i));


            directions = [directions cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'angleSequence'))))];

            tempFR = [tempFR cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'tfSequence'))))];

            spatFR = [spatFR cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'sfSequence'))))];

            interStimStats = cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    Ordered_stims= strsplit(data.VS_ordered{ex},',');
    containsMB = cellfun(@(x) contains(x,'SDG'),Ordered_stims);
    ttlInd = find(containsMB);

     [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,1,"SDG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
     [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,1,"SDG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));

     static_time = cell2mat(stim.VSMetaData.allPropVal(11))*1000; %Static time


     stimInter= mean(stimOn(2:end)-stimOff(1:end-1)); % When dealing with different speeds, save different stim durs, and create M for each speed
     stimDur = mean(-stimOn+stimOff); %When dealing with different speeds, save different stim durs, and create M for each speed
     
     A = [stimOn directions' tempFR' spatFR'];
     [C indexS] = sortrows(A,[2 3 4]);
     %4. Sort directions:
     directimesSorted = C(:,1)';
     B = [stimOff directions' tempFR' spatFR'];
     [Coff indexSo] = sortrows(B,[2 3 4]);
     directimesSortedOff = Coff(:,1)';

     uDir = unique(directions);nDir = length(uDir);
     uTempFR = unique(tempFR);nTempFR = length(uTempFR);
     uSpatFR = unique(spatFR);nSpatFR = length(uSpatFR);
    


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
    preBase = round(stimInter/2);%round(3*interStimStats/4);

    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-stimInter/2)/bin),round((stimDur+stimInter)/bin)); %response matrix
    [nT,nN,nB] = size(Mr);
    trialDivision = nT/(nDir*nSpatFR*nTempFR);

    %figure;imagesc(squeeze(Mr(:,15,:)));colormap(flipud(gray(64)));xline(preBase/bin);xline((preBase+stimDur)/bin);yline(trialDivision:trialDivision:nT-1);

    [MbTotal] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-interStimStats)/bin),round((interStimStats)/bin));% baseline matrix (whole baseline)
    baseline = squeeze(mean(squeeze(mean(MbTotal)),2));
    MrMean = squeeze(mean(Mr));


    stims = stimOn';

    %%%%%%%%%%%%%% Select baseline

    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

    Mb = mean(Mb,3); %mean across time bins

    spkRateBM = mean(Mb); %Calculate baseline spike rate of each neuron.

    %2. Calculate denominator for Z-score calculation
    %spkRateBM = mean(Mb); %total mean.
    epsilon = 0.01;
    %denom = mad(MbC,0)+epsilon; %Calculate baseline variation of each neuron.

    denom = mad(Mb,0)+eps;

    %%%Build rasters

    for Rasters = 1
        rasters = 0;
        if rasters ==1
            preR=500;
            binr = 50;
            [MrRast] =BuildBurstMatrix(goodU,round(p.t/binr),round((directimesSorted-preR)/binr),round((stimDur+preR*2)/binr));
            %response matrix
            [nT,nN,nB] = size(MrRast);

            for u =1:nN
                u=1
                f= figure;
                imagesc(squeeze(MrRast(:,u,:)));colormap(flipud(gray(64)));

            end
        end
    end


    % 1. Statistical test (does the neuron respond to the stimulus in
    % general?)

    %%%% Static
    initialDelay = 100/bin; %100ms of delay
    %%%Given that neurons have different response delays calculate the min
    %%%or max window across the mean trials applying a convolution.
    kernel_size = 300;
    kernel = ones(1,kernel_size) / kernel_size;
    meanTR=squeeze(mean(Mr(:,:,round((stimInter/2))/bin:round((stimInter/2))/bin+static_time), 1));
    % Apply the convolution
    mean_values = zeros(size(meanTR,1),size(meanTR,2)-length(kernel)+1);
    for i = 1:size(meanTR,1)%iterate across neurons
        mean_values(i,:) = conv(meanTR(i,:), kernel, 'valid');
    end

    % Find the maximum mean value and its position
    [staticFR, max_position] = max(mean_values,[],2);
    
    % Adjust the position to account for the 'valid' convolution's offset
    max_position = max_position + floor(kernel_size / 2);

    %%%%
    %mean(Mr(:,:,round((stimInter/2))/bin+initialDelay:round((stimInter/2))/bin+static_time-initialDelay),[3 1]);
    %%To calculate appropiate baseline, find max position in baseline
    %%period (same kernel)
    meanBase=squeeze(mean(Mr(:,:,1:round((stimInter/2))/bin), 1));

    mean_valuesBase = zeros(size(meanBase,1),size(meanBase,2)-length(kernel)+1);
    for i = 1:size(meanTR,1)%iterate across neurons
        mean_valuesBase(i,:) = conv(meanBase(i,:), kernel, 'valid');
    end

    % Find the maximum mean value and its position
    [MaxBase, maxBase_position] = max(mean_valuesBase,[],2);


    %%%
    observed_diffStatic = staticFR - MaxBase;

    %%%% Moving
    movingFR = mean(Mr(:,:,round(preBase/bin)+round(static_time/bin):round(preBase/bin)+round(static_time/bin)+round(stimDur/bin)-round(static_time/bin)-round(preBase/bin)),[3 1]);
    observed_diffMoving = movingFR - spkRateBM;

    rands =1000;
    rand_diffStatic = zeros(rands,nN);
    rand_diffMoving = zeros(rands,nN);

    if ~isfile(sprintf('randValues-It-%d.mat',rands)) || repeatShuff ==1

        for s = 1:rands

            %%%Static: Calculate the specific time window (kernel size) where
            %%%the response is max.
            %response
            M_shuff=Mr(:,:,randperm(size(Mr,3)));
            SHstatic = squeeze(mean(M_shuff(:,:,round((stimInter/2)/bin):round((stimInter/2))/bin+round(static_time/bin)),1));

            mean_valuesSS = zeros(size(SHstatic,1),size(SHstatic,2)-length(kernel)+1);
            for i = 1:size(SHstatic,1)%iterate across neurons
                mean_valuesSS(i,:) = conv(SHstatic(i,:), kernel, 'valid');
            end
            SHstatic = max(mean_valuesSS,[],2);

            %baseline
            Mb_shuffStatic = squeeze(mean(M_shuff(:,:,1:round((stimInter/2))/bin), 1));
            mean_valuesB = zeros(size(Mb_shuffStatic,1),size(Mb_shuffStatic,2)-length(kernel)+1);
            for i = 1:size(Mb_shuffStatic,1)%iterate across neurons
                mean_valuesB(i,:) = conv(Mb_shuffStatic(i,:), kernel, 'valid');
            end
            Mb_shuffStatic = max(mean_valuesB,[],2);

            %%%Moving
            SHmoving = mean(M_shuff(:,:,round(preBase/bin)+round(static_time/bin):round(preBase/bin)+round(static_time/bin)+round(stimDur/bin)-round(static_time/bin)-round(preBase/bin)),[3 1]);
            Mb_shuff = mean(M_shuff(:,:,1:preBase),[3 1]);

            rand_diffStatic(s,:) = SHstatic' - Mb_shuffStatic';
            rand_diffMoving(s,:) = SHmoving - Mb_shuff;
        end

        save(sprintf('SDG-randValues-Static-It-%d',rands),"rand_diffStatic");
        save(sprintf('SDG-randValues-Moving-It-%d',rands),"rand_diffMoving");
    else

        rand_diffStatic = load(sprintf('SDG-randValues-Static-It-%d.mat',rands))
        rand_diffMoving = load(sprintf('SDG-randValues-Moving-It-%d',rands))

    end


    p_valueStat = mean(abs(rand_diffStatic) >= abs(observed_diffStatic)');
    p_valueMov = mean(abs(rand_diffMoving) >= abs(observed_diffMoving));



    %%Sort out neurons that respond less than x% of the trials

    exclusionZone = (trialDivision)/nT;
    sumM = sum(Mr,3);
    exludedN = [];

    inclusionZone = trialDivision*0.7; %Number of trials in the same trial category that are enough as a good unit.  
    
    for u =1:nN

        if sum(sumM(:,u) > 0)<floor(nT*exclusionZone)

            exludedT = ones(1,nT/trialDivision);

            for i=1:trialDivision:nT

                if sum(sumM(i:i+trialDivision-1,u)) <= inclusionZone

                    exludedT(i) = 0;

                end

            end

            if sum(exludedT) >= 1

                exludedN = [exludedN u];

            end

        end

    end

    %%% plot rand diff values as histogram, then real values as lines.
    %%% Then at each side plot a red line where signifficance is
    %%% crossed.

        %1. normalize across neurons random events:

        NstR = normalize(normalize([rand_diffStatic;observed_diffStatic],1),2);

        NmvR = normalize(normalize([rand_diffMoving;observed_diffMoving],1),2);

        %2. normalize observed values across neurons:

        Nst = normalize(observed_diffStatic);
        Nmv = normalize(observed_diffMoving);

        %3. plot random values mean with se bars as histogram

        
        se_valuesST = std(rand_diffStatic, 0, 2) / sqrt(size(rand_diffStatic, 2));
        se_valuesMV = std(rand_diffMoving, 0, 2) / sqrt(size(rand_diffMoving, 2));

        sum((mean(rand_diffStatic,2))>prctile(mean(rand_diffStatic,2), 99.5))/length(mean(rand_diffStatic,2));

        figure;
        h = histogram(mean(NstR(1:end-1,:),2));
        edges = [min(h.BinEdges)-mean(diff(h.BinEdges)), h.BinEdges,max(h.BinEdges)+mean(diff(h.BinEdges))];

        %create plotable observed values
        plotObs = NstR(end,:);
        plotObs(plotObs> max(h.BinEdges)) = max(h.BinEdges)+mean(diff(h.BinEdges))/2;
        plotObs(plotObs<min(h.BinEdges)) = min(h.BinEdges)-mean(diff(h.BinEdges))/2;
        %plot second histogram with observed responses
        hold on; histogram(plotObs,'BinEdges',edges)
        xline(prctile(mean(NstR(1:end-1,:),2), 99.5));xline(prctile(mean(NstR(1:end-1,:),2), 0.5));

        %create plotable observed values

%         pltObs = NstR;
%         plotObs(plotObs>prctile(mean(NstR(1:end-1,:),2), 99.5) ) = 
% 
%         && plotObs>xline(prctile(mean(NstR(1:end-1,:),2), 0.5)
%         xline(NstR(end,:))
        %Produce hist with data such that if the value is above the
        %percentiles, it has a fixed value.
        modObsST = observed_diffStatic;
        %modObsST(observed_diffStatic>) = observed_diffStatic()


        hold on; xline(log10(observed_diffStatic+abs(min(observed_diffStatic))+1)); 
        xline(Nst)
        hold on;
% 
%         % Add error bars for the standard error
%         errorbar(1:length(mean(se_valuesST,2)), mean(se_valuesST,2), se_valuesST, 'k.', 'LineWidth', 1);



        

        %2. normalize 

    % 3. MOVING SECTION

    uDir = unique(directions);
    uDirRad = deg2rad(uDir);
    
    [M3] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted + static_time)/bin),round((stimDur-static_time)/bin)); %For histograms moving


    trialDivision = nT/(length(uDir)); 

    aM = zeros(nN,length(uDir));
    for u =1:nN
        for j = 1:length(uDir)
            aM(u,j) = sum(squeeze(M3(trialDivision*(j-1)+1:trialDivision*j,u,:)),"all"); %Plot cricle histogram.
        end
    end


    OIm = sqrt(sum(aM.*sin(2*uDirRad),2).^2 + sum(aM.*cos(2*uDirRad),2).^2)./sum(aM,2);


    %3.2. SHUFFLED MOVING
    rands =1000;
    OIshM = zeros(nN,rands);
    aSHm = zeros(nN,length(uDir));
    for s = 1:rands

        M_shufM=M3(randperm(size(M3,1)),:,:);

        trialDivision = nT/(length(uDir));

        for u =1:nN
            for j = 1:length(uDir)
                aSHm(u,j) = sum(squeeze(M_shufM(trialDivision*(j-1)+1:trialDivision*j,u,:)),"all"); %Plot cricle histogram.
            end
        end

        OIshM(:,s) = sqrt(sum(aSHm.*sin(2*uDirRad),2).^2 + sum(aSHm.*cos(2*uDirRad),2).^2)./sum(aSHm,2);

    end

    %3.4. MOVING SHUFFLED: PLOT SHUFFLED OI STATISTIC VS LINE WITH REAL

    pmR = (sum(OIshM >= OIm,2) +1)/(rands+1);
    pmL = (sum(OIshM <= OIm) +1)/(rands+1);

    %%Rand
    rand = 1000;

    bin2 =50;

    %%%%Select good responsive neurons:


    


    %%%%%% 

    [Mr2] = BuildBurstMatrix(goodU,round(p.t/bin2),round((directimesSorted-stimInter/2)/bin2),round((stimDur+stimInter)/bin2)); %response matrix
    figure;
    imagesc(squeeze(Mr2(:,1,:)));colormap(flipud(gray(64)));
    %Plot stim start:
    xline(preBase/bin2,'k', LineWidth=1.5)
    %Plot stim end:
    xline(stimDur/bin2+preBase/bin2,'k',LineWidth=1.5)
    ylabel('Trials');xlabel('Time (ms)');
    title(sprintf('U.%d-R.%.3f-B.%.3f',u,max_mean_value(u),spkRateBM(u)));

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





end