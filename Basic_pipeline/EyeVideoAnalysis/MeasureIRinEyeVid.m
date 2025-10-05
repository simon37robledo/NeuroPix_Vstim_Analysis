% Read video
vid = VideoReader('PV132_7_Camera1_20250507-182049.mp4');

% Select a frame number (e.g., frame 100)
frameNum = 800;

% Go to that frame
vid.CurrentTime = (frameNum-1) / vid.FrameRate;
frame = readFrame(vid);

% Show the frame
imshow(frame)
hold on

% Pick two points interactively with the mouse
disp('Click two points in the frame');
[x, y] = ginput(2);   % returns [x1 x2], [y1 y2]

% Plot the selected points
plot(x, y, 'ro-', 'MarkerSize', 10, 'LineWidth', 2)

% Compute Euclidean distance (in pixels)
dist = sqrt( (x(2)-x(1))^2 + (y(2)-y(1))^2 );
title(sprintf('Distance = %.2f pixels', dist))

%%

[angle, lengthEye, midpoint] = eye_line_angle_gui_video('PV132_7_Camera1_20250507-182049.mp4');

deg2rad(angle)