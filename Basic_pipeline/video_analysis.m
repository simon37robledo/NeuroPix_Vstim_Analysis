%% Process video files
j=1;

missingFrames = {};

missedFramesVID ={};

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile,'Format','auto');

%%

for ex = [44]%examplesSDG%[7 8 28]%1:size(data,1)
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

    %2. Extract video dir
    patternIndex = strfind(string(NP.recordingDir), "\catgt");
    endIndex = patternIndex(1)-1;
    vidDir = string(NP.recordingDir);
    vidDir = extractBetween(vidDir,1,endIndex)+"\Original_Video";

    file = dir (vidDir);
    filenames = {file.name};
    vidFrames = filenames(contains(filenames,".csv"));
    vidFrames = vidFrames{contains(vidFrames, "Camera")};

    %vidName = filenames(contains(filenames,".mp4"));


%     figure;
%     plot([1:length(vidTTLs)-1],diff(vidTTLs))

    vidTTLs = Ttrigger{7}; %Camera trigger
    indexTTLs = find(diff(vidTTLs)>1000);

    %vidTTLs(indexTTLs(2))

    cd(vidDir)

    if numel(indexTTLs) > 1 %%Select the triggers within start and stop of video

        vidTTLs = vidTTLs(indexTTLs(1)+1:indexTTLs(2));

    else
        vidTTLs = vidTTLs(indexTTLs(1)+1:end); %If video was not stopped before acquisition ends

    end
    
    VideoTS = readtable(vidDir+filesep+vidFrames); %Real frames saved in video

    timestampsRecF = (VideoTS.timestamp - VideoTS.timestamp(1))*1000; %Zero frames, t-start = 0, and convert to miliseconds

    diffTS = diff(timestampsRecF); %Difference between zeroed frames

    %Plot missing frames
    figure;
    plot([1:size(VideoTS,1)-1],diffTS)

    MissingInd = find(diffTS>15); %find indexes where there is a frame-miss event (could be one or several). 

    %Exclude missing frames in rec time stamps from TTLs

    mFrames= [];

    TTL2FrameI = 1:length(vidTTLs);

    for m = 1:length(MissingInd)

         %Check what is the time period of missed frames within an event and calculate the
         %indexes of the unexisting frames. 

        MissingSection = [MissingInd(m):MissingInd(m)-1+floor(diffTS(MissingInd(m))/10)-1];

        mFrames = [mFrames MissingSection];

    end

    vidTTLsR = vidTTLs;
    vidTTLsR(mFrames) = [];

    TTL2Frame = setdiff(TTL2FrameI,mFrames); %%TTL index number to index number of frame 

    lengthVideoTS = length(VideoTS.timestamp)

    lengthTTLsR = length(vidTTLsR)


    cd(NP.recordingDir)

    save('videoTimeStampsSynced','vidTTLsR')
    save('TTL2Frame','TTL2Frame')

    
%     missedFramesVID{j} =sum(diff(VideoTS.timestamp)>0.015);
% 
%     missingFrames{j} = length(vidTTLsReal)-size(VideoTS,1);
% 
%     missingFramesCorr{j} = length(vidTTLs) - (size(VideoTS,1));

    j=j+1;
    
end

%% Verify with video and full field flashes: %Verified with one video PV35_2

%A. Get times of full field flash


VSordered = strsplit(data.VS_ordered{ex},',');
OBpos = find(VSordered=="OB");

[stimOn stimOff] = NPdiodeExtract(NP,0,0,"NS",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
[stimOn stimOff] = NPdiodeExtract(NP,0,0,"NS",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));

[m indF]=min(abs(stimOn(1)-vidTTLsR));

%triggers TTLs FFF
[m2 indF2]=min(abs(Ttrigger{3}(1)-vidTTLsR));

%B. Get RealTTls frame index and plot frame:
video = VideoReader();

frameNumber = MissingInd; % Example frame number
video.CurrentTime = (indF2 - 1) / video.FrameRate;
frame = readFrame(video);
figure;
imshow(frame)

implay("\\sil3\data\Large_scale_mapping_NP\lizards\PV35\PV35_Experiment_18_8_24\Insertion2\pv35_Insertion2_Camera1_20240819-123015.mp4");


%% Find video sections for moving ball

start = stimOn(1);
stop = stimOff(end);

[~, idxS] = min(abs(timestampsTTLS - start));

[~, idxE] = min(abs(timestampsTTLS - stop));

%% Find movement events

% Define the video file
videoFile = vidDir+filesep+vidName;  % Replace with your actual video file

% Read the video
vidObj = VideoReader(videoFile);

% Get the region of interest (ROI) for cropping
firstFrame = readFrame(vidObj);  % Read the first frame to select the cropping area
imshow(firstFrame);
title('Select the region to analyze');
roi = round(getrect);  % Interactively select the region to crop (ROI)

% Rewind the video to the beginning
vidObj.CurrentTime = 0;

% Initialize variables for motion detection
previousFrame = [];
threshold = 30;  % Adjust this threshold based on the sensitivity of movement detection
motionTimestamps = [];  % Array to store timestamps of detected motion
motionDetected = false;  % Flag to track motion occurrence

% Loop through the frames and detect movement in the ROI
while hasFrame(vidObj)
    % Get the current frame and crop it to the selected ROI
    frame = readFrame(vidObj);  % Read a frame
    croppedFrame = imcrop(frame, roi);  % Crop the frame to the selected ROI
    currentFrame = rgb2gray(croppedFrame);  % Convert to grayscale
    
    currentTime = vidObj.CurrentTime;  % Get the current timestamp in seconds
    
    if ~isempty(previousFrame)
        % Compute the absolute difference between consecutive frames
        frameDifference = abs(double(currentFrame) - double(previousFrame));
        
        % Threshold the difference to detect motion
        motionMask = frameDifference > threshold;
        
        % Count the number of motion pixels
        motionPixels = sum(motionMask(:));
        
        % If motion is detected (based on pixel count threshold), log the timestamp
        if motionPixels > 500  % Adjust this threshold based on the area of movement
            motionTimestamps = [motionTimestamps; currentTime];  % Log timestamp
            motionDetected = true;
        else
            motionDetected = false;
        end
        
        % Display the motion detection frame in the ROI (optional)
        imshow(motionMask);
        title(['Motion Detection in ROI at time: ', num2str(currentTime, '%.2f'), 's']);
        drawnow;
    end
    
    % Update the previous frame
    previousFrame = currentFrame;
end

% Display the motion timestamps
disp('Motion detected at the following timestamps (in seconds):');
disp(motionTimestamps);







