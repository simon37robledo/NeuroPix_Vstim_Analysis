%% Pixels to degreess
%%Dell screen specs

%%Acer active display area 587x330 mm

%function pixels2eyeDegrees()

PixelsX = 1920;

Width = 587;

WidthPix = Width/PixelsX;

angular_size_screen = rad2deg(2 * atan((WidthPix*1500) / (2 * 215)));

angular_distance = rad2deg(atan((WidthPix*750) / 215));


% Known values
x = [1, 2, 3]; % x-coordinates of square corners
y = [4, 5, 6]; % y-coordinates of square corners
z = [7, 8, 9]; % z-coordinates of square corners
L = 10; % Square side length in cm
sphere_center = [0, 0, 0]; % Sphere center coordinates

% Cartesian to Spherical
[azimuth, elevation, radius] = cart2sph(x - sphere_center(1), ...
                                        y - sphere_center(2), ...
                                        z - sphere_center(3));
azimuth = rad2deg(azimuth);
elevation = rad2deg(elevation);

% Angular size
D = radius(1); % Assuming square is at the same distance from the sphere center
angular_size_rad = 2 * atan(L / (2 * D));
angular_size_deg = rad2deg(angular_size_rad);

% Display results
disp('Azimuth (deg):'), disp(azimuth)
disp('Elevation (deg):'), disp(elevation)
disp('Angular size (deg):'), disp(angular_size_deg)




%theta = rad2deg(atan(WidthPix));


