

realTimeFrames = vidTTLsReal/1000-vidTTLsReal(1)/1000;


oddballtimesToVidTimes=(directimesSorted(end-19:end)/1000-vidTTLsReal(1)/1000);

% Load video
video = VideoReader(vidDir+filesep+vidName);

startTimes = oddballtimesToVidTimes+0.2; % Add more start times as needed

% Number of frames to read per start time
numFramesToRead = 5;

numRows = length(startTimes);
k =1;
vidsample = zeros(size(croppedImg,1),size(croppedImg,2),size(croppedImg,3),length(startTimes)*numFramesToRead);
%figure;tiledlayout(length(startTimes),numFramesToRead,"TileSpacing",'none');
% Loop through each start time

video1 = VideoWriter('output_video.mp4', 'MPEG-4');
video1.FrameRate = 30; % Adjust frame rate as needed
open(video1);
for i = 1:length(startTimes)
    startTime = startTimes(i);
    startFrame = round(startTime * video.FrameRate);

    % Read and plot two frames per start time
    for j = 0:numFramesToRead-1
        frameIndex = startFrame + 20;
        %subplot(numRows, numFramesToRead, (i-1)*numFramesToRead + j + 1);
        if frameIndex <= video.NumFrames
            frame = read(video, frameIndex);  % Read frame       
            % Display frame
            croppedImg = imcrop(frame, rect);
            vidsample(:,:,:,k) = croppedImg;
            writeVideo(video1, croppedImg);
            k = k+1;
            %title(['Start Time: ' num2str(startTime) ' s, Frame: ' num2str(frameIndex)]);
        end
    end
end

close(video1)
implay(vidsample)

size(vidsample)

%%
frameNumber = 100; % Example frame number
video.CurrentTime = (frameNumber - 1) / video.FrameRate;
frame = readFrame(video);

% Display the frame and open cropping GUI
imshow(frame);
[croppedFrame rect] = imcrop;
