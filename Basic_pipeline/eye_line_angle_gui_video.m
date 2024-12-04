function angle = eye_line_angle_gui_video(video_file)
    % Input: video_file (string) - Path to the video file
    % Output: angle (double) - Final selected angle of the line

    % Check if the video file exists
    if ~isfile(video_file)
        error('The specified video file does not exist.');
    end

    % Create the GUI figure
    fig = figure('Name', 'Eye Angle Calculator (Video Frame)', ...
                 'NumberTitle', 'off', ...
                 'MenuBar', 'none', ...
                 'ToolBar', 'none', ...
                 'Position', [100, 100, 600, 600]);
    
    % Load the video file
    vid = VideoReader(video_file);
    
    % Extract the first frame (or any specific frame)
    frame_number = 1; % Modify to get a specific frame
    vid.CurrentTime = (frame_number - 1) / vid.FrameRate; % Set to the desired frame
    frame = readFrame(vid);
    
    % Display the frame in the GUI
    ax = axes('Parent', fig, 'Position', [0.1, 0.3, 0.8, 0.6]);
    imshow(frame, 'Parent', ax);
    title(ax, 'Draw and Adjust the Line, Then Select the Angle', 'FontSize', 12);

    % Initialize the angle variable
    angle = NaN;
    is_selected = false; % To check if the angle is finalized
    
    % Allow the user to draw one line
    h = imline(ax);
    if isempty(h)
        disp('No line drawn. Exiting.');
        close(fig);
        return;
    end

    % Callback to update the angle dynamically when the line moves
    addNewPositionCallback(h, @(pos) update_angle(pos));
    
    % Add a text box to display the angle
    angle_display = uicontrol('Style', 'text', ...
                              'String', 'Angle: N/A', ...
                              'Position', [150, 30, 300, 30], ...
                              'FontSize', 12, ...
                              'HorizontalAlignment', 'left');
    
    % Add a button to finalize the angle
    uicontrol('Style', 'pushbutton', ...
              'String', 'Select Angle', ...
              'Position', [480, 30, 100, 30], ...
              'Callback', @select_angle_callback);
    
    % Wait for the user to finalize the angle
    uiwait(fig);

    % Nested function to calculate and update the angle
    function update_angle(pos)
        % Get the positions of the line endpoints
        x1 = pos(1, 1);
        y1 = pos(1, 2);
        x2 = pos(2, 1);
        y2 = pos(2, 2);

        % Calculate the angle (in degrees) of the line with respect to the x-axis
        dx = x2 - x1;
        dy = y2 - y1;
        angle = atan2d(dy, dx); % Angle in degrees

        % Update the angle display
        angle_display.String = sprintf('Angle: %.2f°', angle);
    end

    % Callback to finalize the angle
    function select_angle_callback(~, ~)
        % Set the flag to true and close the GUI
        is_selected = true;
        uiresume(fig);
        close(fig);
    end
end
