function rotated_output = rotateMatrix(input_matrix, pivot_point, angle_degrees, interpolation_method, show_visualization)
% ROTATEMATRIX Rotates a 2D matrix or image around a specified pivot point
%
%   ROTATED_OUTPUT = ROTATEMATRIX(INPUT_MATRIX, PIVOT_POINT, ANGLE_DEGREES) 
%   rotates the INPUT_MATRIX around PIVOT_POINT by ANGLE_DEGREES using 
%   'linear' interpolation
%
%   ROTATED_OUTPUT = ROTATEMATRIX(INPUT_MATRIX, PIVOT_POINT, ANGLE_DEGREES, INTERPOLATION_METHOD)
%   allows specifying the interpolation method ('nearest', 'linear', or 'cubic')
%
%   ROTATED_OUTPUT = ROTATEMATRIX(INPUT_MATRIX, PIVOT_POINT, ANGLE_DEGREES, INTERPOLATION_METHOD, SHOW_VISUALIZATION)
%   controls whether to display a visualization of the rotation (true/false)
%
%   Inputs:
%       INPUT_MATRIX     - 2D matrix or image to rotate
%       PIVOT_POINT      - [x, y] coordinates of rotation center
%       ANGLE_DEGREES    - Rotation angle in degrees (positive = counterclockwise)
%       INTERPOLATION_METHOD - Optional: 'nearest', 'linear', or 'cubic' (default: 'linear')
%       SHOW_VISUALIZATION  - Optional: true or false (default: false)
%
%   Output:
%       ROTATED_OUTPUT   - Rotated matrix with same dimensions as input

    % Handle optional arguments
    if nargin < 4
        interpolation_method = 'linear';
    end
    
    if nargin < 5
        show_visualization = false;
    end

    % Validate inputs
    validateattributes(input_matrix, {'numeric'}, {});  % Just check it's numeric, no dimension restriction
    validateattributes(pivot_point, {'numeric'}, {'vector', 'numel', 2});
    validateattributes(angle_degrees, {'numeric'}, {'scalar'});

    % Check if input is probably an image (has 3 dimensions for RGB)
    is_image = ndims(input_matrix) == 3 || (ismatrix(input_matrix) && max(input_matrix(:)) <= 1);

    % Create the affine2d transformation
    tform = affine2d();
    %counter clock-wise rotation:
%     tform.T = [cos(angle_degrees) -sin(angle_degrees) 0;
%                sin(angle_degrees) cos(angle_degrees) 0;
%                0 0 1];

% For clockwise rotation with positive angles
tform.T = [cos(angle_degrees) sin(angle_degrees) 0;
          -sin(angle_degrees) cos(angle_degrees) 0;
           0 0 1];
    
    % Set the rotation center (pivot point)
    tform.T(3,1:2) = pivot_point - pivot_point*tform.T(1:2,1:2);
    
    % Apply transformation with specified interpolation method
    rotated_output = imwarp(input_matrix, tform, 'OutputView', imref2d(size(input_matrix)), ...
                            'InterpolationMethod', interpolation_method);
    
    % Show visualization if requested
    if show_visualization
        figure;
        
        % First subplot - original matrix/image
        subplot(1, 2, 1);
        if is_image
            imshow(input_matrix);
            title('Original Image');
        else
            imagesc(input_matrix);
            colorbar;
            title('Original Matrix');
            colormap('jet');
        end
        
        hold on;
        plot(pivot_point(1), pivot_point(2), 'ro', 'MarkerSize', 10, 'LineWidth', 2);
        hold off;
        
        % Second subplot - rotated matrix/image
        subplot(1, 2, 2);
        if is_image
            imshow(rotated_output);
            title('Rotated Image');
        else
            imagesc(rotated_output);
            colorbar;
            title('Rotated Matrix');
            colormap('jet');
        end
        
        sgtitle(['Rotation: ' num2str(angle_degrees) '° using ' interpolation_method ' interpolation']);
    end
end