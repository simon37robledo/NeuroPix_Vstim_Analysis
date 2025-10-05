
ImagePath = "\\sil3\data\Large_scale_mapping_NP\lizards\PV132\PV132_Experiment_4_5_25\Images_setup";

cd(ImagePath);


NP.recordingDir
files = dir;
names = string({files.name});
IRfile = names(contains(names,"IRmirror"));
img = imread(ImagePath+filesep+IRfile);
figure;imagesc(img);

Xc = 636;
Yc = 1056;

Xm1 = 2298;
Ym1 = 1247;
Xm2 = 2475;
Ym2 = 1077;

Xe = 2429;
Ye = 1040;


Xc1 = -16.36;
Yc = -3.379;

Xm1 = 4.542;
Ym1 = -4.078;

Xm2 = 6.102;
Ym2 = -2.355;

Xe = 5.077;
Ye = -2.012;
%%
%PV132 
% imagesPaths = {'\\sil3\data\Large_scale_mapping_NP\lizards\PV132\PV132_Experiment_4_5_25\Images_setup\4_5_25_In1-2\IRmirror_angle_pic.jpg',...
%     '\\sil3\data\Large_scale_mapping_NP\lizards\PV132\PV132_Experiment_4_5_25\Images_setup\5_5_25_In3-4\IRmirror_angle_in3-4-Cam.jpg',...
%     '\\sil3\data\Large_scale_mapping_NP\lizards\PV132\PV132_Experiment_4_5_25\Images_setup\6_5_25_In5-6\IRmirror_angle_In5-6.jpg',...
%      '\\sil3\data\Large_scale_mapping_NP\lizards\PV132\PV132_Experiment_4_5_25\Images_setup\7_5_25_In7-9\IR_mirror_image_In_7-9.jpg'};
%ins = {[1,2],[3,4],[5,6],[7,9]};
% PV140
% imagesPaths = {'\\sil3\data\Large_scale_mapping_NP\lizards\PV140\PV140_Experiment_4_6_25\Images_setup\4_6_25_In1-3\IR_mirror_image_1-3.jpg',...
%     '\\sil3\data\Large_scale_mapping_NP\lizards\PV140\PV140_Experiment_4_6_25\Images_setup\5_6_25_In4-6\IR_mirror_image_4-6.jpg',...
%     '\\sil3\data\Large_scale_mapping_NP\lizards\PV140\PV140_Experiment_4_6_25\Images_setup\9_6_25_In7-11\IR_mirror_Image_In7-11.jpg'};
%ins = {[1,3],[4,6],[7,11]};

%PV51
imagesPaths = {'\\sil3\data\Large_scale_mapping_NP\lizards\PV51\PV51_Experiment_20_8_25\CraniotomyAndsetupPics\UpPicIns-1-3.jpg',...
    '\\sil3\data\Large_scale_mapping_NP\lizards\PV51\PV51_Experiment_20_8_25\CraniotomyAndsetupPics\UpPictureIns-4-5.jpg'};
ins = {[1,3],[4,5]};

%cd('\\sil3\data\Large_scale_mapping_NP\lizards\PV140\PV140_Experiment_4_6_25\Images_setup');
cd('\\sil3\data\Large_scale_mapping_NP\lizards\PV51\PV51_Experiment_20_8_25\CraniotomyAndsetupPics')
%coorTable = table();


for i = 2%2:numel(imagesPaths)

    iPath = imagesPaths{i};

    fc = coordinateSelectorGUI(iPath);
    fc.Animal =  repmat({'PV140'}, height(fc), 1);
    fc.Session = repmat({ins{i}}, height(fc), 1);
    coorTable = [coorTable; fc]; 

end

writetable(coorTable,'coordinates.csv');

%%
j = imread(iPath);

figure;imagesc(j);

size(j,1)


%%

coorTable = readtable('coordinates.csv');

kx = 0;%[0:20]-10;

%Row 1 PV132In1-2 | R2 PV132In3-4 | R2 PV132In5-6 | R2 PV132In7-9

FusionCoorsX = [-34.069,7.35,11.23,9.25;-16.36,4.54,6.102,5.077;-17.103,4.747,6.03,4.576;-23.87,10.173,13.572,11.265];%Xc, Xm1, Xm2, Xe
FusionCoorsY = [-3.858,-8.42,-4.746,-3.6;-3.379,-4.078,-2.355,-2.012;2.123,-0.03,1.3,1.98;2.133,-0.78,2.844,3.074];%Yc, Ym1, Ym2, Ye

Fusion =0;

for i = 1:height(coorTable)

    XmoV = [];
    aV = [];
    Xm2V =[];
    Xm1V =[];
    XciDV =[];

    fig = figure;

    for k = 1:length(kx)

        if Fusion

            Xc = FusionCoorsX(i,1)+kx(k);%coorTable.Camera_X(i)/10;%
            Yc = FusionCoorsY(i,1);%(size(j,1)-coorTable.Camera_Y(i))/-10;%

            Xm1 =  FusionCoorsX(i,2)+kx(k);%(coorTable.Mirror1_X(i))/10;%
            Ym1 = FusionCoorsY(i,2);%(size(j,1)-coorTable.Mirror1_Y(i))/-10;%

            Xm2 =  FusionCoorsX(i,3)+kx(k);%coorTable.Mirror2_X(i)/10;%
            Ym2 = FusionCoorsY(i,3);%(size(j,1)-coorTable.Mirror2_Y(i))/-10;%

            Xe =  FusionCoorsX(i,4)+kx(k); %coorTable.Eye_X(i)/10;%
            Ye =  FusionCoorsY(i,4); %(size(j,1)-coorTable.Eye_Y(i))/-10;%
        else
            Xc = coorTable.Camera_X(i)+kx(k);%coorTable.Camera_X(i)/10;%
            Yc = coorTable.Camera_Y(i);%(size(j,1)-coorTable.Camera_Y(i))/-10;%

            Xm1 = coorTable.Mirror1_X(i)+kx(k);%(coorTable.Mirror1_X(i))/10;%
            Ym1 = coorTable.Mirror1_Y(i);%(size(j,1)-coorTable.Mirror1_Y(i))/-10;%

            Xm2 = coorTable.Mirror2_X(i)+kx(k);%coorTable.Mirror2_X(i)/10;%
            Ym2 = coorTable.Mirror2_Y(i);%(size(j,1)-coorTable.Mirror2_Y(i))/-10;%

            Xe =  coorTable.Eye_X(i)+kx(k); %coorTable.Eye_X(i)/10;%
            Ye =  coorTable.Eye_Y(i); %(size(j,1)-coorTable.Eye_Y(i))/-10;%

        end

        %%%Camera
        f = sqrt((Ym1-Yc).^2+(Xm1-Xc).^2);

        b = sqrt((Ye-Yc).^2+(Xe-Xc).^2);

        m = (Ym2 - Ym1)/(Xm2 - Xm1);
        n = Ym2-m*Xm2;

        syms Xmo a

        d=sqrt((Xmo-Xc).^2+(m*Xmo+n - Yc).^2);
        c=sqrt((Xmo-Xe).^2+(m*Xmo+n - Ye).^2);
        e=sqrt((Xmo-Xm1).^2 + (m*Xmo+n - Ym1).^2);

        Y = vpasolve([b^2 == d^2 + c^2 - 2*d*c*cos(2*a),...
            f.^2 == d^2 + e^2 - 2*e*d*cos(pi/2-a)],[Xmo,a],[Xm1+Xm1*0.01 Xm2-Xm2*0.01; 0.2 1.5]);


        syms Xci

        Xmo=double(Y.Xmo);
        Ymo=double(Y.Xmo)*m+n;
        d=sqrt((Xmo-Xc).^2+(m*Xmo+n - Yc).^2);


        m1 = (Ye - Ymo)/(Xe - Xmo);
        n1 = Ye-m1*Xe;

        Y1 = vpasolve((Xci-Xmo)^2+(Xci*m1+n1-Ymo)^2 == d^2,Xci);

        XciD = double(Y1);
        YciD = double(m1).*XciD+double(n1);

        %     coorTable.NewCam_X(i) = XciD;
        %     coorTable.NewCam_Y(i) = YciD;

        XmoV = [XmoV Xmo];
        aV = [aV double(Y.a)];
        Xm1V = [Xm1V Xm1];
        Xm2V = [Xm2V Xm2];
        XciDV = [XciDV XciD];

    end

    [my,ind] = min([abs(YciD(1)-Xm1), abs(YciD(2)-Xm1)]);
    coorTable.newCamX(i) = XciD(ind);
    coorTable.newCamY(i) = YciD(ind);
    coorTable.Cam_eye_angle(i) = double(Y.a);


    subplot(3,1,1)
    scatter(kx,XmoV);
    title(string(i))
    hold on
    scatter(kx,Xm1V,'red','filled');
    scatter(kx,Xm2V,'yellow','filled')
    xlabel('Constant added to X coor')
    ylabel('X mirror-cam intersection')

    subplot(3,1,2)
    scatter(kx,aV);
    ylabel('Angle')
    xlabel('Constant added to X coor')

    subplot(3,1,3)
    scatter(kx,XciDV(1,:));
    hold on
    scatter(kx,XciDV(2,:));
    ylabel('New camera position')
    xlabel('Constant added to X coor')
      
    figure;
    scatter([Xc;Xm1;Xm2;Xe],[Yc;Ym1;Ym2;Ye],100,[1:4],"filled");hold on;
    text([Xc;Xm1;Xm2;Xe],[Yc;Ym1;Ym2;Ye],{'Cam','Mirror1','Mirror2','Eye'});
    plot([0 10],[0*m+n,10*m+n]);
    plot(double(Xmo),double(Ymo),'*k')
    title(string(i))
    plot(XciD,YciD,'ok')
    axis equal;


end
%% Convert into cm and save

coorTableTransf = coorTable;

for i = 1:height(coorTable)
    distBoard(i) = sqrt(sum(([coorTable.Edge1_board1_X(i) coorTable.Edge1_board1_Y(i) ] - [coorTable.Edge2_board_X(i)  coorTable.Edge2_board_Y(i) ]).^2));
end

distReal = 30;
factor = distReal./distBoard;

for i = 1:height(coorTable)
coorTableTransf(i,[1:14 18:20]) = array2table(coorTable{i,[1:14 18:20]}*factor(i));
end

%%Use eye center as 0,0

coorTableTransf.newCamX = coorTableTransf.newCamX - coorTableTransf.Eye_X;
coorTableTransf.newCamY = coorTableTransf.newCamY - coorTableTransf.Eye_Y;

writetable(coorTableTransf,'coordinatesTransf.csv');

%% Check cam angle
% Define points as vectors
camN = [XciD(ind), YciD(ind)];  % Replace with actual coordinates
Mirror = [Xmo, Ymo];
Eye = [Xe,Ye];

% Compute vectors
u = Mirror - camN;  % vector AB
v = Eye - Mirror;  % vector BC

% Compute the dot product and norms
dot_uv = dot(u, v);
norm_u = norm(u);
norm_v = norm(v);

% Compute angle in radians
angle_rad = acos(dot_uv / (norm_u * norm_v));

% Convert to degrees
angle_deg = rad2deg(angle_rad);

% Display the result
fprintf('Angle between AB and BC is %.2f degrees\n', angle_deg);


%% Check solutions

a = double(Y.a);
Xmo = double(Y.Xmo);


f = sqrt((Ym1-Yc).^2+(Xm1-Xc).^2);

b = sqrt((Ye-Yc).^2+(Xe-Xc).^2);

m = (Ym2 - Ym1)/(Xm2 - Xm1);
n = Ym2-m*Xm2;

d=sqrt((Xmo-Xc).^2+(m*Xmo+n - Yc).^2);
c=sqrt((Xmo-Xe).^2+(m*Xmo+n - Ye).^2);
e=sqrt((Xmo-Xm1).^2 + (m*Xmo+n - Ym1).^2);

Ymo=Xmo*m+n;

m1 = (Ye - Ymo)/(Xe - Xmo);
n1 = Ye-m1*Xe;

Xci = double(Y1);
Yci = double(m1).*Xci+double(n1);


b^2 == d^2 + c^2 - 2*d*c*cosd(2*a)

f.^2 == d^2 + e^2 - 2*e*d*cosd(90-a)

(Xci(2)-Xmo)^2+(Xci(2)*m1+n1-Ymo)^2 == d^2





%%
writetable(coorTable,'coordinates.csv');

%%
figure;
scatter([Xc;Xm1;Xm2;Xe],[Yc;Ym1;Ym2;Ye],100,[1:4],"filled");hold on;
text([Xc;Xm1;Xm2;Xe],[Yc;Ym1;Ym2;Ye],{'Cam','Mirror1','Mirror2','Eye'});
plot([0 10],[0*m+n,10*m+n]);
plot(Xmo,Ymo,'*k')
plot(Xci,Yci,'ok')
axis equal;

%% 
function coordsTable = coordinateSelectorGUI(imgPath)
    % Validate image path
    if nargin < 1 || ~isfile(imgPath)
        error('Please provide a valid image path.');
    end

    % Load image
    img = imread(imgPath);
    
    %Rotate image if portrait (height > width)
    if size(img, 1) > size(img, 2)
        img = permute(img, [2 1 3]);  % Rotate 90° counterclockwise without interpolation
    end

    % Output placeholder
    coordsTable = [];

    % Create UI figure
    fig = uifigure('Name', 'Select Coordinates', 'Position', [100 100 900 700]);
    ax = uiaxes(fig, 'Position', [50 100 800 550]);
    title(ax, 'Click to select 4 points: Camera, Eye, Mirror1, Mirror2');
    
    imshow(img, 'Parent', ax);
    
    hold(ax, 'on');

    % Instructions
    uilabel(fig, 'Position', [50 660 500 20], ...
        'Text', 'Click on image to select: Camera, Eye, Mirror Edge 1, Mirror Edge 2, IR center, Edge1 board1, Edge2 board');

    % Initialize data
    pointNames =  {'Camera', 'Eye', 'Mirror1', 'Mirror2', 'IR_light','Edge1_board1', 'Edge2_board'};
    draggablePoints = gobjects(1,numel(pointNames));
    coords = nan(numel(pointNames),2);
    pointCount = 0;

    % Accept button
    btn = uibutton(fig, 'Text', 'Accept', ...
        'Position', [750, 660, 100, 30], ...
        'Enable', 'off', ...
        'ButtonPushedFcn', @(btn,event) onAccept());

    % Click listener
    fig.WindowButtonDownFcn = @onClick;

    % Pause execution until Accept is clicked
    uiwait(fig);

    % Click to place points
    function onClick(~, ~)
        if pointCount >= numel(pointNames)
            return;
        end
        cp = ax.CurrentPoint;
        x = cp(1,1);
        y = cp(1,2);
        if x < 0 || y < 0 || x > size(img,2) || y > size(img,1)
            return;
        end
        pointCount = pointCount + 1;
        flippedY = size(img, 1) - y;
        %coords(pointCount,:) = [x, flippedY];
        coords(pointCount,:) = [x,y];
        draggablePoints(pointCount) = drawpoint(ax, ...
            'Position', [x y], ...
            'Color', 'r', ...
            'InteractionsAllowed', 'translate', ...
            'Label', pointNames{pointCount}, ...
            'LabelVisible', 'on');
        if pointCount == numel(pointNames)
            btn.Enable = 'on';
        end
    end

% Accept: collect data, close GUI, return table
    function onAccept()
        for i = 1:numel(pointNames)
            pos = draggablePoints(i).Position;
            %coords(i,:) = [pos(1), size(img, 1) - pos(2)];  % Flip Y
            coords(i,:) = [pos(1), pos(2)];  
        end
        % Flatten and relabel as single-row table
        flatCoords = reshape(coords', 1, []);  % [X1 Y1 X2 Y2 ...]

        % Build column names: Camera_X, Camera_Y, etc.
        varNames = cell(1, numel(pointNames)*2);
        for i = 1:numel(pointNames)
            varNames{2*i - 1} = [pointNames{i}, '_X'];
            varNames{2*i}     = [pointNames{i}, '_Y'];
        end

        coordsTable = array2table(flatCoords, 'VariableNames', varNames);

        uiresume(fig);
        close(fig);
    end

end

