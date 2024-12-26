% Parameters for the raster plot and arrow plotting
nT = 80; % Number of trials
nB = 50; % Number of time points
direcN = 4; % Number of direction types (up, left, down, right)


uDir = unique(directions);
% up, left, down, right

figure;
% Raster plot
subplot('Position', [0.05, 0.1, 0.7, 0.8]); % Adjust to leave space for arrows
imagesc(rand(nT, nB)); % Random data for raster plot
colormap gray;
xlabel('Time');
ylabel('Trials');
title('Raster Plot');
%set(gca, 'YDir', 'reverse'); % Reverse y-axis

dirStart = C(1,2);
offStart = C(1,3);
for t = 1:nT*mergeTrials
    if dirStart ~= C(t,2)
        yline(t/mergeTrials+0.5,'r',LineWidth=3);
        dirStart = C(t,2);
    end
    if offStart ~= C(t,3)
        yline(t/mergeTrials+0.5,'b',LineWidth=2);
        offStart = C(t,3);
    end

end

%%%%%%
% Get the current axis handle for the raster plot
% Get the current axis handle for the raster plot
ax = gca;

setX = 10;
% Parameters for arrow plotting
arrowXStart = nB - setX; % Position the arrows at the last 10 columns of the raster plot
arrowLength = 2;      % Default length of the arrow
headSize = 1.5;       % Size of the arrowhead

uDir =sort(uDir,'descend');
% Plot the arrows for each trial type (every 20 trials per type)
for i = 1:direcN
    % Get the direction from uDir array
    direction = uDir(i);
    halftrials = nT/direcN/2;
    
    % Modify arrow length for up/down arrows to make the shaft longer
    if direction == 0  % Up or Down
        arrowLength = halftrials/2;ys = halftrials-1;
    elseif direction == pi 
        arrowLength = halftrials/2;ys = -halftrials/2;
    else
        arrowLength = 2;ys=0;  % Default length for left/right arrows
    end
    
    % The starting position of the arrow for each trial type
    if direction == pi/2  % Right
        startX = nB - (setX-1);
        ys = 1;% Start at the last 10 columns for right arrows
        direction = direction +pi/2;
    elseif direction == 3*pi/2  % Left
        startX = nB -setX-1;
        direction = direction +pi/2;
        ys = -(3*halftrials)/4;% Start at the last 10 columns minus arrow length for left arrows
    else  % Up or Down
        startX = nB - setX; 
         direction = direction +3*pi/2;
        % Use the same starting X position for up and down
    end

   

    startY = (i - 1) * (nT/direcN) + (nT/direcN) / 2 - ys;% Increase the length for up/down arrows

    % For vertical arrows (up or down), position them in the middle of the trials
     % Centered in the trial block
    
    % Call the plotArrow function to plot the arrows
    plotArrow(ax, direction, arrowLength, startX, startY, headSize);
end

% Function to plot arrows
function plotArrow(ax, direction, arrowLength, arrowXStart, arrowYStart, headSize)
    % Normalize the arrow's direction
    dx = arrowLength * cos(direction);  % Calculate dx based on direction (cosine for x-direction)
    dy = -arrowLength * sin(direction);  % Calculate dy based on direction (negative sine for y-direction)
    
    % Get the axis limits to normalize the coordinates
    xLimits = ax.XLim;
    yLimits = ax.YLim;

    % Convert the starting position from data to normalized figure coordinates
    normXStart = (arrowXStart - xLimits(1)) / (xLimits(2) - xLimits(1));
    normYStart = (arrowYStart - yLimits(1)) / (yLimits(2) - yLimits(1));

    % Calculate the end coordinates of the arrow in normalized figure space
    normXEnd = (arrowXStart + dx - xLimits(1)) / (xLimits(2) - xLimits(1));
    normYEnd = (arrowYStart + dy - yLimits(1)) / (yLimits(2) - yLimits(1));
    
    % Plot the arrow using normalized coordinates (to be between 0 and 1)
    annotation('arrow', ...
        [normXStart, normXEnd], ...  % X positions in normalized figure coordinates
        [normYStart, normYEnd], ...  % Y positions in normalized figure coordinates
        'LineWidth', 2, ...
        'HeadWidth', headSize * 2, ...
        'HeadLength', headSize * 2, ...
        'HeadStyle', 'cback2', ...
        'Color', 'k');  % Black color for the arrow
end



