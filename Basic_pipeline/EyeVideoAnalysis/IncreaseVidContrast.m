% Input video

patternIndex = strfind(string(NP.recordingDir), "\catgt");
endIndex = patternIndex(1)-1;
ellipseDir= string(NP.recordingDir);
ellipseDir = extractBetween(ellipseDir,1,endIndex);
file = dir (ellipseDir);
filenames = {file.name};

input_video = filenames(contains(filenames,".mp4"));
output_video = sprintf('EnhancedContrast-%s',input_video{2});

% Read video
cd(ellipseDir)
video_reader = VideoReader(input_video{2}); %Select cropped
video_writer = VideoWriter(output_video, 'MPEG-4');
open(video_writer);


% Define cropping rectangle (e.g., [x, y, width, height])
%crop_rect = [100, 50, 640, 480]; % Modify as needed

while hasFrame(video_reader)
    frame = readFrame(video_reader);
    
    % Convert to grayscale (optional, depends on your contrast adjustment needs)
    gray_frame = rgb2gray(frame);
    
    % Enhance contrast
    enhanced_frame = imadjust(gray_frame,[0 0.3],[]);
       
    % Convert to RGB (if needed for saving video)
    output_frame = repmat(enhanced_frame, [1, 1, 3]); 
    
    % Write frame to output video
    writeVideo(video_writer, output_frame);
end

% Close video writer
close(video_writer);

disp('Video processing complete.');