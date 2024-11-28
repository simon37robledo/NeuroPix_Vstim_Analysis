

%function EyePositionAnalysis(NP,divisions,varargin)%plots,newRun)
divisions =11;
newRun =1;
plots=1;
% Mandatory Inputs:
%   NP - Neuropixels class
%   divisions - Number of divisions for the grid. 
% Optional Inputs (varargin):
%   varargin{1} (plots)  - Plot ellipse centers on top of closest to mean
%   ellipse center frame.
%   varargin{2} (stimType) - n by three matrix where the first column is
%   onset time, 2nd is offset time and third is stim type index. If false,
%   colors are assigned depending of position on grid.
%   varargin{3} (newRun) - Rewrite previous run with same parameters
% Default values for optional parameters
% plots = false;
% stimType = false;
% newRun = false;

% % Assign optional inputs if provided
% if nargin > 2
%     plots = varargin{1};
% end
% if nargin > 3
%     stimType = varargin{2}; %n by three matrix where the first column is onset time, 2nd is offset time and third is stim type index
% end
% if nargin > 4
%     newRun = varargin{3};% Rewrite previous run with same parameters
% end

if mod(divisions,2)~=1
    disp('Divisions can not be even in order to obtain a central quadrant')
    return;
end

if ~exist(sprintf('timeSnipsNoMov-%d-quadrants.mat',divisions), 'file') || newRun
patternIndex = strfind(string(NP.recordingDir), "\catgt");
endIndex = patternIndex(1)-1;
ellipseDir= string(NP.recordingDir);
ellipseDir = extractBetween(ellipseDir,1,endIndex);
file = dir (ellipseDir);
filenames = {file.name};

ellipseName = filenames(contains(filenames,"Eye_ellipse_thr"));
goodFrames = filenames(contains(filenames, "good_frames_thr"));
videoFile =  filenames(contains(filenames,".mp4"));

%ellipsePath = "\\sil3\data\Large_scale_mapping_NP\lizards\PV35\PV35_Experiment_18_8_24\Insertion2\PV35_2_Eye_ellipse";

Data = readtable(ellipseDir+filesep+ellipseName); %Check which rows are Nans with the elipse code
video = VideoReader(ellipseDir+filesep+videoFile);
goodFrames = readtable(ellipseDir+filesep+goodFrames);


%%Calculate the mean of the center of the eye and then calculate a
%%rectangle for cropping whose center square contains the mean center of
%%the eye

meanCenter = [mean(Data.center_x), mean(Data.center_y)];

%%% Find frame number that is closest to mean:
[mn indEC] = min(abs(Data.center_x-meanCenter(1))+abs(Data.center_y-meanCenter(2)));
%indEC = 1;
frameNumber = goodFrames.Values(indEC); % Example frame number

video.CurrentTime = (frameNumber - 1) / video.FrameRate;
frame = readFrame(video);

% Display the frame and open cropping GUI
%imshow(frame);
% [croppedFrame rect] = imcrop;    %rect = [xmin ymin width height]
% 
% [mx ind] = max([rect(3),rect(4)]);

%Find max square possible for the center mean

maxSizeX = size(frame,2)/2;
maxSizeY = size(frame,1)/2;

if meanCenter(1) > maxSizeX
    maxSquareX =  maxSizeX*2-meanCenter(1);
else
    maxSquareX =  meanCenter(1);
end

if meanCenter(2) > maxSizeY
    maxSquareY =  maxSizeY*2-meanCenter(2);
else
    maxSquareY =  meanCenter(2);
end 

maxSquareSize = 2*round((min(maxSquareX,maxSquareY)));

rectNew = [meanCenter(1)-maxSquareSize/2 meanCenter(2)-maxSquareSize/2 maxSquareSize maxSquareSize];

croppedFrameC = imcrop(frame, rectNew);
% 
% figure;
% imshow(croppedFrameC);


%%%Get indexes of frames in which the center of the eye that falls in different subdivisions of the screen. 

%EyeLoc =cell(2,divisions*divisions);
grid = linspace(0, maxSquareSize, divisions + 1);
j=1;
indFrames= zeros(1, length(Data.center_x),'single');
for x = 1:divisions
    for y =1:divisions
        xIndx = find((Data.center_x > rectNew(1)+grid(x)) & (Data.center_x < (rectNew(1) + grid(x) + mean(diff(grid)))));
        yIndx = find((Data.center_y > rectNew(2)+grid(y)) & (Data.center_y < (rectNew(2) + grid(y) + mean(diff(grid)))));
        %EyeLoc{1,j} =  intersect(xIndx, yIndx);
        indFrames(intersect(xIndx, yIndx)) = zeros(1,length(intersect(xIndx, yIndx)))+j;
        j =j+1;
    end
end

%%%Ephys sync times to frames
cd(NP.recordingDir)
frameTimes = load('videoTimeStampsSynced.mat').vidTTLsR;
frameTimes = frameTimes(goodFrames.Values);

%%%%Within a square id, find time intervals

indexesChange=find(diff(indFrames)~=0)+1;
timeSnips = zeros(3,length(indexesChange)+1);

%%Asign start and end timeSnip
timeSnips(:,1) = [0;frameTimes(indexesChange(1));indFrames(1)];
timeSnips(:,end) = [frameTimes(indexesChange(end));frameTimes(end);indFrames(end)];

for i = 1:length(indexesChange)-1
    timeSnips(1:2,i+1) = [frameTimes(indexesChange(i));frameTimes(indexesChange(i+1))];
    timeSnips(3,i+1) = indFrames(indexesChange(i));
end


save(sprintf('timeSnipsNoMov-%d-quadrants',divisions),'timeSnips')

if plots ==1
%%%Plot data
% Example data
centers = [Data.center_x, Data.center_y]; % Replace with your X, Y data
[X, Y] = meshgrid(grid, grid);

figure;
imagesc(croppedFrameC)
% xlim([0 max(grid)])
% ylim([0 max(grid)])
axis equal

hold on;
plot(X, Y, 'w', X', Y', 'w');
% xlim([0 max(grid)]);
% ylim([0 max(grid)]);% Plot grid lines
xticks(grid)
yticks(grid)
axis equal;

if length(stimType)> 1

    posXPerTrial = [];
    posYPerTrial = [];
    indexType = [];
    trialIndex = [];
    %get frameTimes that are between the beginning and end of stim
    for i =1:length(stimType)
        posXPerTrial = [posXPerTrial Data.center_x(frameTimes>=stimType(i,1) & frameTimes<stimType(i,2))'];
        posYPerTrial = [posYPerTrial Data.center_y(frameTimes>=stimType(i,1) & frameTimes<stimType(i,2))'];
        indexType = [indexType zeros(1,length(frameTimes(frameTimes>=stimType(i,1) & frameTimes<stimType(i,2))))+stimType(i,3)];
        trialIndex = [trialIndex zeros(1,length(frameTimes(frameTimes>=stimType(i,1) & frameTimes<stimType(i,2))))+i];
    end

    % Border colors for each group
    borderColorsV = zeros(length(indexType),3);
    borderColorsV(indexType'==max(indexType),:) = repmat([0, 0, 0.5],sum(indexType'==max(indexType)),1);
    borderColorsV(indexType'==min(indexType),:) = repmat([0.5, 0, 0],sum(indexType'==min(indexType)),1);

    colorVector =  (stimType(:,4) - min(stimType(:,4))) / (max(stimType(:,4)) - min(stimType(:,4)));
    colorVector = colorVector(trialIndex);

    % Define the lightened versions of colors
    C_light = min(borderColorsV + 0.8, 1); % Add lightness, ensuring values stay <= 1

    % Interpolate between dark and light colors using B
    CombinedColors = borderColorsV.*colorVector + C_light.*(1 - colorVector); % Blend based on B
    

    scatter(posXPerTrial-rectNew(1), posYPerTrial-rectNew(2),25,CombinedColors,'filled') %'CData', colorVector', 'filled');

    hold on;
    for i = unique(indexType)
        scatter(x(group==i), y(group==i), 50, 'MarkerFaceColor', colors(i,:), ...
            'MarkerEdgeColor', 'none', 'MarkerFaceAlpha', 0.5);
    end
    hold off;

    scatter(posXPerTrial-rectNew(1), posYPerTrial-rectNew(2),borderColors);
    title(sprintf('Eye Movements %s',NP.recordingName));
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    axis equal;

else
    gscatter(Data.center_x-rectNew(1), Data.center_y-rectNew(2),indFrames);
    title(sprintf('Eye Movements %s',NP.recordingName));
    xlabel('X Coordinate');
    ylabel('Y Coordinate');
    axis equal;
end

scatter(meanCenter(1)-rectNew(1),meanCenter(2)-rectNew(2),25, 'filled','MarkerFaceColor','r')
axis equal;
%scatter(Data.center_x(indEC)-rectNew(1),Data.center_y(indEC)-rectNew(2),25, 'filled','MarkerFaceColor','w')
axis equal;
% Parameters of the ellipse
a = Data.width(indEC); % Semi-major axis
b = Data.height(indEC); % Semi-minor axis
h = Data.center_x(indEC)-rectNew(1); % X-coordinate of the center
k = Data.center_y(indEC)-rectNew(2); % Y-coordinate of the center
phi =  Data.phi(indEC); % Rotation angle in radians (45 degrees)
theta = 0:0.01:2*pi; % Angle range for the ellipse

% Parametric equations for the ellipse with rotation
x = h + a * cos(theta) * cos(phi) - b * sin(theta) * sin(phi);
y = k + a * cos(theta) * sin(phi) + b * sin(theta) * cos(phi);
% Plot the ellipse
plot(x, y, 'r', 'LineWidth', 2);
axis equal;
legend('off')
end

else

    fprintf('Analysis done! saved in: timeSnipsNoMov-%d-quadrants.mat',divisions)

end

%end








