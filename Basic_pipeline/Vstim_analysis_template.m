%% Vstim analysis function template


cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%SA5_1,PV103_1,PV27_1

for ex =[1 8 28]%1:size(data,1)

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


 %2. Extract moving stim statistics
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

     [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,1,"StaticDriftingGratings",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));

end