
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