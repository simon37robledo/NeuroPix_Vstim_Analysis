%% Process video files
j=1;

missingFrames = {};

missedFramesVID ={};

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



    [Ttrigger,chNumberT]=NP.getTrigger();


    endIndex = patternIndex(1)-1;
    vidDir = string(NP.recordingDir);
    vidDir = extractBetween(vidDir,1,endIndex);

    file = dir (vidDir);
    filenames = {file.name};

    file = dir (vidDir);
    filenames = {file.name};
    vidFrames = filenames(contains(filenames,".csv"));

    vidTTLs = Ttrigger{7};

    figure;
    plot([1:length(vidTTLs)-1],diff(vidTTLs))

    indexTTLs = find(diff(vidTTLs)>1000);

    if numel(indexTTLs) > 1

        vidTTLsReal = vidTTLs(indexTTLs(1)+1:indexTTLs(2));

    else
        vidTTLsReal = vidTTLs(indexTTLs(1)+1:end);

    end

   
    VideoTS = readtable(vidDir+filesep+vidFrames);

     %Check timing of missing frames:
    figure;
    plot([1:size(VideoTS,1)-1],diff(VideoTS.timestamp))

    diffTS = diff(VideoTS.timestamp);
    
    numMissedFrames = diffTS(diffTS>0.015)*100;

    timestampsTTLS= (vidTTLsReal-vidTTLsReal(1));

    timestampsRecF = (VideoTS.timestamp - VideoTS.timestamp(1))*1000;

    %Exclude missing frames in rec time stamps from TTLs

    diffTS = diff(timestampsRecF);

    MissingInd = find(diffTS>15);

    mFrames= [];

    for m = 1:length(MissingInd)

        MissingSection = [MissingInd(m):MissingInd(m)-1+floor(diffTS(MissingInd(m))/10)-1];

        mFrames = [mFrames MissingSection];
    end


    timestampsTTLS(mFrames) = [];

    cd(NP.recordingDir)

    save('videoTimeStamps','timestampsTTLS')


    missedFramesVID{j} =sum(diff(VideoTS.timestamp)>0.015);

    missingFrames{j} = length(vidTTLsReal)-size(VideoTS,1);

    missingFramesCorr{j} = length(timestampsTTLS) - (size(VideoTS,1));

    j=j+1;
end




