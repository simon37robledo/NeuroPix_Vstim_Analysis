path_brainO = '\\132.66.45.127\data\Large_scale_mapping_NP\Model_brain\PV_Brain_original_DVR_.stl';

Image2d =  '\\132.66.45.127\data\Large_scale_mapping_NP\Imaging\29-9-22-PVslices\Snap-58.jpeg'; 

% Load the 3D model
brain = stlread(path_brainO);

% Load the 2D image
texture = imread(Image2d);

% Create a figure
figure;

% Create a surface mesh
h = patch("Faces",brain.ConnectivityList,"Vertices",brain.Points,'FaceLighting',   'gouraud', ...
    'EdgeColor',       'none', ...
    'FaceColor',       [0.8 0.8 1.0], ...
    'AmbientStrength', 0.15, 'FaceAlpha', 0.3);

% Map the 2D image as a texture
set(h, 'FaceColor', 'texturemap', 'CData', texture, 'FaceAlpha', 1);

% Set axis limits and labels
axis equal; % Ensure aspect ratio is maintained
xlabel('X');
ylabel('Y');
zlabel('Z');

% Optionally, save the rendered 3D object with the applied texture as an image
print('output_image.png', '-dpng', '-r300'); % Save as PNG with 300 DPI