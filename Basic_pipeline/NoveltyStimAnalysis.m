cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

controlNS = 1;

%%
for ex = [40 41 42 43]%examplesSDG%[7 8 28]%1:size(data,1)
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
    rectFiles = filenames(contains(filenames,"rectGrid"));

    
        if isempty(rectFiles)
            %disp()
            w= sprintf('No rect grid files where found in %s. Skipping into next experiment.',NP.recordingName);
            warning(w)
            continue
        end


    directions = [];
    offsets = [];
    sizes = [];
    speeds = [];
    orientations = [];

    j =1;

    VSordered = strsplit(data.VS_ordered{ex},',');
    RGpos = find(VSordered=="RG");
    OBpos = find(VSordered=="OB");
    OBCpos = find(VSordered=="OBC");

    [orderVS orderVSIndex] = sort([RGpos OBpos OBCpos]);

     selecFiles = rectFiles(orderVSIndex(2));

     if controlNS

         selecFiles = rectFiles(orderVSIndex(3));
     end

   

    positions = [];
  %  NSfiles = 
    if size(rectFiles) ~= [0 0]

        for i = selecFiles %Extract stim properties
            ball= load(stimDir+"\"+string(i));

            positions = [positions cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'pos'))))];

            interStimStats = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end


    %3. Load Triggers (diode)
    ttlInd = OBpos;

    if controlNS
        ttlInd = OBCpos;
    end

    if controlNS

        [stimOn stimOff ] = NPdiodeExtract(NP,0,"NSC",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
        [stimOn stimOff] = NPdiodeExtract(NP,0,"NSC",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A

    else
        [stimOn stimOff ] = NPdiodeExtract(NP,0,"NS",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
        [stimOn stimOff] = NPdiodeExtract(NP,0,"NS",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A
    end

    stimInter= mean(stimOn(2:end)-stimOff(1:end-1)); % When dealing with different speeds, save different stim durs, and create M for each speed
    stimDur = mean(-stimOn+stimOff); %When dealing with different speeds, save different stim durs, and create M for each speed

    A = [stimOn positions'];

    [C indexS] = sortrows(A,[2]);

    directimesSorted = C(:,1)';

%     LFP = NP.getData(300,round(stimOn-stimInter/2),round(stimDur+stimInter));
% 
%     LFP = LFP(:,indexS',:);
% 
% 
%     %size(LFP,3)/(NP.samplingFrequency/1000)
% 
%     figure;
%     imagesc(squeeze(LFP));caxis([-3000 10000])
% 
%     figure;
%     imagesc(squeeze(LFP(100,:,:)));caxis([-200 1200])


    %5. Load data to select good units
    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');
    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");
    GoodU_orDepth = cluster_info.depth(cluster_info.group=="good");

    %If new tic matrix needs to be used, move oldTIc files and run convertPhySorting2tIc
    cd(NP.recordingDir)

    %Load KS results
    p = NP.convertPhySorting2tIc(NP.recordingDir);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');
    goodSPKS = p.nSpks(label == 'good');
    goodAmp = p.neuronAmp(label == 'good');

    %7. Load raster matrices
    bin = 20;
    preBase = round(stimInter/2);%round(3*interStimStats/4);

    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((stimDur+preBase*2)/bin)); %response matrix
    [nT,nN,nB] = size(Mr);

    % Preallocate normalized matrix
    nMr = zeros(nT, nN, nB);

    uPos = unique(positions);


    if ~controlNS
        for u = 1:nN
            f = figure('Visible', 'off');

            imagesc(squeeze(Mr(:,u,:)));colormap(flipud(gray(64)));caxis([0 1])
            title(sprintf('NoveltyStim|%s|UnitN-%d|UnitPhy-%d',strrep(NP.recordingName,'_','-'),u,GoodU_or(u)))

            xline(preBase/bin, '-g', LineWidth=1.5);
            xline(preBase/bin+stimDur/bin, '-b', LineWidth=1.5);
            yline(nT-min([sum(positions==uPos(1)) sum(positions==uPos(2))]),'-','Odd position',LineWidth=1.5,LabelHorizontalAlignment="left")
            ylabel('Trials')
            xlabel('time (ms)')
            xticks([0:round((preBase/2)/bin):size(Mr,3)])
            xticklabels([0:round((round(preBase/2/10))*10):size(Mr,3)*bin])
            set(f,'Color','w')
            cd(NP.recordingDir+"\Figs")
            f.Position =[680   224   560   754];
            print(f, sprintf('NoveltyStim-%s-Unit-%d.png',NP.recordingName,u),'-dpng');
            close

        end
    else

        posUnames = num2cell(string(uPos));
        for u = 1:nN
            f = figure('Visible', 'on');

            imagesc(squeeze(Mr(:,u,:)));colormap(flipud(gray(64)));caxis([0 1])
            title(sprintf('NoveltyStimControl|%s|UnitN-%d|UnitPhy-%d',strrep(NP.recordingName,'_','-'),u,GoodU_or(u)))

            xline(preBase/bin, '-g', LineWidth=1.5);
            xline(preBase/bin+stimDur/bin, '-b', LineWidth=1.5);
            yline([nT/length(uPos):nT/length(uPos):nT],'-',posUnames,LabelHorizontalAlignment="right")
            ylabel('Trials')
            xlabel('time (ms)')
            xticks([0:round((preBase/2)/bin):size(Mr,3)])
            xticklabels([0:round((round(preBase/2/10))*10):size(Mr,3)*bin])
            set(f,'Color','w')
            cd(NP.recordingDir+"\Figs")
            f.Position =[680   224   560   754];
            print(f, sprintf('NoveltyStimControl-%s-Unit-%d.png',NP.recordingName,u),'-dpng');
            close

        end
        
    end


    %%%Control


end