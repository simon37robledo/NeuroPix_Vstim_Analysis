
video = VideoReader('C:\Users\MarkS9\Desktop\HeadFixedEyeTrack2-Simon-2024-10-02\EyeTracking3-Simon-2024-10-03\videos\pv35_Insertion3_Camera1_20240819-152541DLC_Resnet50_EyeTracking3Oct3shuffle1_snapshot_200_p60_labeled.mp4');

outputVideo = VideoWriter('PV35_3_output_with_ellipse.mp4', 'MPEG-4');
open(outputVideo);

Ellipse_data = readtable('filename.csv');

while hasFrame(video)
    frame = readFrame(video);

    % Plot frame and ellipse
    imshow(frame); hold on;
    ellipse_x = 100; % x-coordinate of ellipse center
    ellipse_y = 100; % y-coordinate of ellipse center
    ellipse_width = 50; % width of ellipse
    ellipse_height = 30; % height of ellipse
    angle = 45; % rotation angle in degrees
    
    % Plot ellipse
    ellipseHandle = ellipse(ellipse_width, ellipse_height, deg2rad(angle), ellipse_x, ellipse_y, 'r');
    
    % Capture the frame with the ellipse
    frameWithEllipse = getframe(gca);
    writeVideo(outputVideo, frameWithEllipse.cdata); % write frame if saving video

    hold off;
end