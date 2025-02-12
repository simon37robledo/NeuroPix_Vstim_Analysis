function [angle, length, midpoint] = eye_line_angle_gui_video(video_file)
    % Input: video_file (string) - Path to the video file
    % Outputs:
    %   angle (double) - Final angle of the line
    %   length (double) - Length of the line
    %   midpoint (1x2 vector) - [x, y] coordinates of the line's midpoint

    % Check if the video file exists
    if ~isfile(video_file)
        error('The specified video file does not exist.');
    end

    % Create the GUI figure
    fig = figure('Name', 'Eye Line Calculator (Video Frame)', ...
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
    title(ax, 'Draw and Adjust the Line, Then Select It', 'FontSize', 12);

    % Initialize outputs
    angle = NaN;
    length = NaN;
    midpoint = [NaN, NaN];
    is_selected = false; % To check if the line is finalized
    
    % Allow the user to draw one line
    h = imline(ax);
    if isempty(h)
        disp('No line drawn. Exiting.');
        close(fig);
        return;
    end

    % Callback to update the angle, length, and midpoint dynamically when the line moves
    addNewPositionCallback(h, @(pos) update_line_details(pos));
    
    % Add a text box to display the angle, length, and midpoint
    details_display = uicontrol('Style', 'text', ...
                                 'String', 'Angle: N/A | Length: N/A | Midpoint: N/A', ...
                                 'Position', [50, 30, 500, 30], ...
                                 'FontSize', 10, ...
                                 'HorizontalAlignment', 'left');
    
    % Add a button to finalize the line
    uicontrol('Style', 'pushbutton', ...
              'String', 'Select Line', ...
              'Position', [480, 30, 100, 30], ...
              'Callback', @select_line_callback);
    
    % Wait for the user to finalize the line
    uiwait(fig);

    % Nested function to calculate and update line details
    function update_line_details(pos)
        % Get the positions of the line endpoints
        x1 = pos(1, 1);
        y1 = pos(1, 2);
        x2 = pos(2, 1);
        y2 = pos(2, 2);

        % Calculate the angle (in degrees) of the line with respect to the x-axis
        dx = x2 - x1;
        dy = y2 - y1;
        angle = atan2d(dy, dx); % Angle in degrees
        
        % Calculate the length of the line
        length = sqrt(dx^2 + dy^2);
        
        % Calculate the midpoint of the line
        midpoint = [(x1 + x2) / 2, (y1 + y2) / 2];

        % Update the details display
        details_display.String = sprintf('Angle: %.2f° | Length: %.2f | Midpoint: [%.2f, %.2f]', ...
                                         angle, length, midpoint(1), midpoint(2));
    end

    % Callback to finalize the line
    function select_line_callback(~, ~)
        % Set the flag to true and close the GUI
        is_selected = true;
        uiresume(fig);
        close(fig);
    end
end