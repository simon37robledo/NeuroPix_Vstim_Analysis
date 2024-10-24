%% Process video files
j=1;

missingFrames = {};

missedFramesVID ={};

for ex = 41%examplesSDG%[7 8 28]%1:size(data,1)
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
    vidName = filenames(contains(filenames,".mp4"));

    vidTTLs = Ttrigger{7};

    figure;
    plot([1:length(vidTTLs)-1],diff(vidTTLs))

    indexTTLs = find(diff(vidTTLs)>1000);

    if numel(indexTTLs) > 1 %%Select the triggers within start and stop of video

        vidTTLsReal = vidTTLs(indexTTLs(1)+1:indexTTLs(2));

    else
        vidTTLsReal = vidTTLs(indexTTLs(1)+1:end);

    end

   
    VideoTS = readtable(vidDir+filesep+vidFrames); %actually saved frames in video

    %Check timing of missing frames:
    figure;
    plot([1:size(VideoTS,1)-1],diff(VideoTS.timestamp))

    diffTS = diff(VideoTS.timestamp);
    
    numMissedFrames = diffTS(diffTS>0.015)*100;

    %%% Modify digital triggers to reflect real frames
    length(diffTS)+round(numMissedFrames);

    length(vidTTLsReal)

    timestampsTTLS= vidTTLsReal;%(vidTTLsReal-vidTTLsReal(1));

    timestampsRecF = (VideoTS.timestamp - VideoTS.timestamp(1))*1000;

    %Exclude missing frames in rec time stamps from TTLs

    diffTS = diff(timestampsRecF);

    MissingInd = find(diffTS>15);

    mFrames= [];

    TTL2FrameI = 1:length(timestampsTTLS);

    for m = 1:length(MissingInd)

        MissingSection = [MissingInd(m):MissingInd(m)-1+floor(diffTS(MissingInd(m))/10)-1];

        mFrames = [mFrames MissingSection];

    end


    TTL2Frame = setdiff(TTL2FrameI,mFrames);


    timestampsTTLS(mFrames) = [];

    cd(NP.recordingDir)

    save('videoTimeStamps','timestampsTTLS')
    save('TTL2Frame','TTL2Frame')

    

    missedFramesVID{j} =sum(diff(VideoTS.timestamp)>0.015);

    missingFrames{j} = length(vidTTLsReal)-size(VideoTS,1);

    missingFramesCorr{j} = length(timestampsTTLS) - (size(VideoTS,1));

    j=j+1;
end

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







