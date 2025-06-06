%% Process video files
% 1. Find good and synced frames from the digital triggers that correspond to real captured frames in reptilearn
% 2. Filter Jitter of DLC results in case there is. 

j=1;

missingFrames = {};

missedFramesVID ={};

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile,'Format','auto');

%%

%% Process video files
j=1;

missingFrames = {};

missedFramesVID ={};

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile,'Format','auto');

removeVibrations =0;


%%

for ex = [57:72]%examplesSDG%[7 8 28]%1:size(data,1):66
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class

    NP = loadNPclassFromTable(ex);

    [Ttrigger,chNumberT]=NP.getTrigger();
    
    vidDir = data.Eye_video_dir{ex};


    file = dir (vidDir);
    filenames = {file.name};
    %Original video path:
    vidFrames = filenames{contains(filenames,".csv") & contains(filenames,"Camera") & ~contains(filenames,"snapshot") & contains(filenames,"_"+string(data.Insertion(ex))+"_")};


    vidTTLs = Ttrigger{7}; %Camera trigger
    indexTTLs = find(diff(vidTTLs)>1000);

    %vidTTLs(indexTTLs(2))

    cd(vidDir)

    if numel(indexTTLs) > 1 %%Select the triggers within start and stop of video

        vidTTLs = vidTTLs(indexTTLs(1)+1:indexTTLs(2));

    else
        vidTTLs = vidTTLs(indexTTLs(1)+1:end); %If video was not stopped before acquisition ends

    end

    VideoTS = readtable(string(vidDir)+filesep+string(vidFrames)); %Real frames saved in video

    timestampsRecF = (VideoTS.timestamp - VideoTS.timestamp(1))*1000; %Zero frames, t-start = 0, and convert to miliseconds

    diffTS = diff(timestampsRecF); %Difference between zeroed frames


    msPerFarme = mean(diffTS);
    frameRate = round(1000/msPerFarme);

    if ex == 65 || ex == 66 ||  ex == 67
        removeVibrations = 1;
    else
        removeVibrations = 0;
    end


    if removeVibrations

        filename = filenames{contains(filenames,".csv") & contains(filenames,"snapshot") & contains(filenames,"_"+string(data.Insertion(ex))+"_")};
        T =  readtable(filename, ...
            'ReadVariableNames', true, ...
            'NumHeaderLines', 2);  % Skip first 3 lines to get to the actual data;

        fs = frameRate;  % Frames per second (adjust based on your video)
        f0 = 49;              % Frequency to remove (Hz)
        Q = 10;               % Quality factor (higher = narrower notch)

        % Design notch filter
        wo = f0 / (fs/2);     % Normalized frequency
        bw = wo / Q;
        [b, a] = iirnotch(wo, bw);

        % Create a copy to store filtered data
        T_filtered = T;

        % Get all variable (column) names
        varNames = T.Properties.VariableNames;

        % Loop through all columns
        for i = 1:length(varNames)
            colName = varNames{i};

            % Check if column is x or y coordinate
            if contains(colName, 'x') || contains(colName, 'y')
                Tempdata = T{:, i};


                % Apply zero-phase lowpass filter (ignore NaNs)
                nanIdx = isnan(Tempdata);
                dataClean = Tempdata;
                dataClean(nanIdx) = interp1(find(~nanIdx), Tempdata(~nanIdx), find(nanIdx), 'linear', 'extrap');
                filtered = filtfilt(b, a, dataClean);

                % Restore NaNs after filtering
                filtered(nanIdx) = NaN;

                % Save into filtered table
                T_filtered{:, i} = filtered;
            end
        end

%         %%%Test
%         figure;plot(2000:2500,filtered(2000:2500)')
%         title('Notch, Q = 5')
%         figure;plot(2000:2500,T.y(2000:2500)')

        % Path to original file and output file
        originalFile = filename;
        outputFile   = sprintf('%s-FilteredJitter.csv',extractBefore(filename, ".csv"));

        % --- Step 1: Read original headers ---
        fid = fopen(originalFile);
        header1 = fgetl(fid);  % Line 1: scorer
        header2 = fgetl(fid);  % Line 2: bodypart
        header3 = fgetl(fid);  % Line 3: coords (x/y/likelihood)
        fclose(fid);

        % --- Step 2: Write headers to output file ---
        fid_out = fopen(outputFile, 'w');
        fprintf(fid_out, '%s\n', header1);
        fprintf(fid_out, '%s\n', header2);
        fprintf(fid_out, '%s\n', header3);
        fclose(fid_out);

        % --- Step 3: Append the filtered data ---
        % Write numeric data (no variable names, just values)
        writetable(T_filtered, outputFile, ...
            'WriteVariableNames', false, ...
            'WriteMode', 'append');


    end

    MissingInd = find(diffTS>(msPerFarme+5));

    mFrames= [];

    TTL2FrameI = 1:length(vidTTLs);


    for m = 1:length(MissingInd)

        %Check what is the time period of missed frames within an event and calculate the
        %indexes of the unexisting frames.

        MissingSection = [MissingInd(m):MissingInd(m)-1+floor(diffTS(MissingInd(m))/msPerFarme)-1];

        mFrames = [mFrames MissingSection];

    end

    vidTTLsR = vidTTLs;
    vidTTLsR(mFrames) = [];

    if isempty(mFrames)
        mFrames = 0;
    end

   fprintf('Missed frames = %d',length(mFrames))

    TTL2Frame = setdiff(TTL2FrameI,mFrames); %%TTL index number to index number of frame

    lengthVideoTS = length(VideoTS.timestamp);

    lengthTTLsR = length(vidTTLsR);
    %
%     figure;plot(vidTTLsR(mFrames(1)-50:mFrames(1)-50+200-1)-vidTTLsR(1),ones(1,200),'*')
%     test = vidTTLsR(mFrames(1)-50:mFrames(1)-50+200-1)-vidTTLsR(1);
%     xline(timestampsRecF(mFrames(1)-50:mFrames(1)-50+200-1));


    framesInsec = sum(timestampsRecF >= timestampsRecF(1) & timestampsRecF<=timestampsRecF(1)+1000);

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







