cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%%
for ex = 40%examplesSDG%[7 8 28]%1:size(data,1)
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

    [stimOn stimOff ] = NPdiodeExtract(NP,0,"RG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff] = NPdiodeExtract(NP,0,"RG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A


    stimInter= mean(stimOn(2:end)-stimOff(1:end-1)); % When dealing with different speeds, save different stim durs, and create M for each speed
    stimDur = mean(-stimOn+stimOff); %When dealing with different speeds, save different stim durs, and create M for each speed

    A = [stimOn positions'];

    [C indexS] = sortrows(A,[2]);

    LFP = NP.getData(1:384,round(stimOn-stimInter/2),round(stimDur+stimInter));

    LFP = LFP(:,indexS',:);

    %size(LFP,3)/(NP.samplingFrequency/1000)

    figure;
    imagesc(squeeze(mean(LFP,1)))

    figure;
    imagesc(squeeze(LFP(100,:,:)));caxis([-200 1200])


    %5. Load data to select good units
    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');
    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");
    GoodU_orDepth = cluster_info.depth(cluster_info.group=="good");

    %If new tic matrix needs to be used, move oldTIc files and run convertPhySorting2tIc
    cd(NP.recordingDir)




end