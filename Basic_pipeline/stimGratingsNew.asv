%% stimGratings

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

newTIC = 0;
%%
%SA5_1,PV103_1,PV27_1

for ex =[2 8 28]%1:size(data,1)

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


    % 1. Statistical test (does the neuron respond to the stimulus in
    % general?)

    %%%% Static
    initialDelay = 100/bin; %100ms of delay
    staticFR = mean(Mr(:,:,round((stimInter/2))/bin+initialDelay:round((stimInter/2))/bin+static_time-initialDelay)/(static_time-initialDelay),[3 1]);
    observed_diffStatic = staticFR - spkRateBM;

    %%%% Moving
    movingFR = mean(Mr(:,:,round(preBase/bin)+round(static_time/bin):round(preBase/bin)+round(static_time/bin)+round(stimDur/bin)-round(static_time/bin))-round(preBase/bin),[3 1]);
    observed_diffMoving = movingFR - spkRateBM;



    rands =1000;
    rand_diffStatic = zeros(s,nN);
    rand_diffMoving = zeros(s,nN);

        for s = 1:rands

            M_shuff=Mr(:,:,randperm(size(Mr,3)));
            SHstatic = mean(M_shuff(:,:,round((stimInter/2)/bin)+initialDelay:round((stimInter/2))/bin+round(static_time/bin)-initialDelay),[3 1]);

            SHmoving = mean(M_shuff(:,:,round(preBase/bin)+round(static_time/bin):round(preBase/bin)+round(static_time/bin)+round(stimDur/bin)-round(static_time/bin)...
                -round(preBase/bin)),[3 1]);

            Mb_shuff = mean(M_shuff(:,:,1:preBase),[3 1]);
            
            rand_diffStatic(s,:) = SHstatic - Mb_shuff;
            rand_diffMoving(s,:) = SHmoving - Mb_shuff;
        end
        p_valueStat = mean(abs(rand_diffStatic) >= abs(observed_diffStatic));
        p_valueMov = mean(abs(rand_diffMoving) >= abs(observed_diffMoving));

        %%% plot rand diff values as histogram, then real values as lines.
        %%% Then at each side plot a red line where signifficance is
        %%% crossed. 


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