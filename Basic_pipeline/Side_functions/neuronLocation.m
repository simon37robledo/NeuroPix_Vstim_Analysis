% ex=14;
% data_ex = data(ex,:);
% 
% plotN = 1;
% 
% plotIns =1;

function [Coor, Fig] = neuronLocation(NP,data_ex,goodU,plotN,plotIns,colors,Fig,newCoor)
%%%NP = neuropixels class
%%%goodU = TIC matrix
%%%PlotN = plot neurons (logical)
%%%PlotIns = plot insertions (logical)
%%%Fig = figure of 3D brain on top of whcih to plot neurons or insertions.


if ~exist(string(NP.recordingDir)+filesep+string(NP.recordingName)+"_g0_tcat.imec0.ap_kilosortChanMap.mat")

    g =SGLXMetaToCoords(string(NP.recordingDir)+filesep+string(NP.recordingName)+"_g0_tcat.imec0.ap.meta");

end

chanMap = load(string(NP.recordingName)+"_g0_tcat.imec0.ap_kilosortChanMap.mat");

path_Brain = '\\sil1\data\Large_scale_mapping_NP\Model_brain\PV_103_imaged_inclined_45TA.stl';

brain = stlread(path_Brain);

if (plotN || plotIns) && isempty(Fig) 
    close all

    Fig = figure;
    axes('Parent', gcf);
    hold on;

    patch("Faces",brain.ConnectivityList,"Vertices",brain.Points,'FaceLighting',   'gouraud', ...
        'EdgeColor',       'none', ...
        'FaceColor',       [0.8 0.8 1.0], ...
        'AmbientStrength', 0.15, 'FaceAlpha', 0.3);

    camlight('headlight');
    material('dull');
    % Fix the axes scaling, and set a nice view angle
    axis('image');
    sideview = [90 0];
    topview = [-90 90];
    frontview = [0 0];
    view(frontview)

    Ylims = ylim;
    Xlims = xlim;
    Zlims = zlim;

    ymax = 80;
    zmin =-65;

    ylim([Ylims(1)+15 ymax])
    zlim([zmin Zlims(2)])

    % yticks([5:5:ymax-30])
    zticks([zmin:5:round(Zlims(2))])
    xticks([round(Xlims(1)):5:round(Xlims(2))])

    grid on
    xlabel('Y (mm)')
    ylabel('X (mm)')
    zlabel('Z (mm)')
    xticklabels(0:0.5:Xlims(2)/10);
    yticklabels(0.5:0.5:ymax/10);
    zticklabels(round((zmin - Zlims(2)))/10:0.5:0);
    sizeDot = 50;
    sizeText = 13;
    offsetTextX = -1.3;
    offsetTextY = 0;
    offsetTextZ = 1;

elseif plotN || plotIns

    figure(Fig)
    sideview = [90 0];
    topview = [-90 90];
    frontview = [0 0];

end


%%%Load coordinates from excel data

if ~isfile(sprintf('Insertion_coordinates-%s.mat',NP.recordingName)) || newCoor

    if string(cell2mat(data_ex.Probe)) == "NP2.0"

        if length(cell2mat(data_ex.Matlab3Dcoor))>3 %%Check if the insertion has matlab coordinates

            parts = split(cell2mat(data_ex.Matlab3Dcoor), ';');

            Xs = split(parts{1},',');
            x1 = str2double(Xs{1});
            x2 = str2double(Xs{2});

            Ys =  split(parts{2},',');
            y1 = str2double(Ys{1});
            y2 = str2double(Ys{2});

        else
            x1 = 0;
            x2 = 0;

            y1 = 0;
            y2 = y1 - 750/100;
        end

        % Vectors for AB and AC
        AB = [0, y2 - y1];
        AC = [x2 - x1, y2 - y1];

        % Calculate the dot product and magnitudes of AB and AC
        dotProduct = dot(AB, AC);
        magAB = norm(AB);
        magAC = norm(AC);

        % Calculate the angle in radians and convert to degrees
        angleRad = acos(dotProduct / (magAB * magAC));

        j=1;

        coors = zeros(length(unique(chanMap.kcoords)),3);
        for i = unique(chanMap.kcoords)'

            %If shank id is:
            if i==1
                coors(j,:) = [x1,y1,surfaceZ(x1,y1, brain)];
            end

            if i==2
                coors(j,:)= [x1 - 750/3/100 * cos(angleRad), y1 - 750/3/100 * sin(angleRad), surfaceZ(x1 - 750/3/100 * cos(angleRad), y1 - 750/3/100 * sin(angleRad),brain)];
            end

            if i==3
                coors(j,:) = [x1 - 750/(3/2)/100 * cos(angleRad), y1 - 750/(3/2)/100 * sin(angleRad), surfaceZ(x1 - 750/(3/2)/100 * cos(angleRad), y1 - 750/(3/2)/100 * sin(angleRad),brain)];

            end

            if i==4

                coors(j,:) = [x1 - 750/100 * cos(angleRad), y1 - 750/100 * sin(angleRad), surfaceZ(x1 - 750/100 * cos(angleRad), y1 - 750/100 * sin(angleRad),brain)];

            end

            j=j+1;

            %                 coorS1 = [x1,y1,surfaceZ(x1,y1, brain)];
            %                 coorS4 = [x1 - 750/100 * cos(angleRad), y1 - 750/100 * sin(angleRad), surfaceZ(x1 - 750/100 * cos(angleRad), y1 - 750/100 * sin(angleRad),brain)];
            %                 coorS2 = [x1 - 750/3/100 * cos(angleRad), y1 - 750/3/100 * sin(angleRad), surfaceZ(x1 - 750/3/100 * cos(angleRad), y1 - 750/3/100 * sin(angleRad),brain)];
            %                 coorS3 = [x1 - 750/(3/2)/100 * cos(angleRad), y1 - 750/(3/2)/100 * sin(angleRad), surfaceZ(x1 - 750/(3/2)/100 * cos(angleRad), y1 - 750/(3/2)/100 * sin(angleRad),brain)];
        end
        x = coors(:,1)';
        y =  coors(:,2)';
        z = coors(:,3)';

        j =1;
        shankIds = unique(chanMap.kcoords);

        uX = zeros(2,size(goodU,2));
        uY = zeros(2,size(goodU,2));
        uZ = zeros(2,size(goodU,2));

        for u = 1:size(goodU,2)


            %Main unit channel position

            %change channel 0 to 1

            ch = goodU(1,u);

            uX(1,u) = x(shankIds == chanMap.kcoords(ch)); %STL coor
            uX(2,u) = x(shankIds == chanMap.kcoords(ch))*100; %Standard coor in uM

            uYvert = y(shankIds == chanMap.kcoords(ch));
            uZvert = data_ex.coor_Z - chanMap.ycoords(ch); %Vertical depth in uM

            uZ(1,u) = max(z)-sin(deg2rad(data_ex.Angle))*(uZvert/100); %depth of unit along vertical axis in STL coor
            uZ(2,u) = sin(deg2rad(data_ex.Angle))*(uZvert); %depth of unit along vertical axis zeroed in uM

            uY(1,u) = uYvert+cos(deg2rad(data_ex.Angle))*(uZvert/100); %Y  of unit along in STL coor
            uY(2,u) = (uYvert*100)+cos(deg2rad(data_ex.Angle))*(uZvert); %Y  of unit along in standard coor in uM

            if plotN
                if ~isempty(Fig) 
                    hold on
                end
                try
                scatter3(uX(1,u),uY(1,u),uZ(1,u),25,colors(u,:),'filled')
                catch
                    colors(u,:)
                end
                view(sideview)
            end

            j=j+1;


        end

        Coor = {uX,uY,uZ;... %unit coordinates
            [x;x],... %X shank coordinates STL
            [y;y+cos(deg2rad(data_ex.Angle))*(data_ex.coor_Z/100)],... %Y shank coordinates
            [repmat(max(z),[1 2]);repmat(max(z),[1 2])-sin(deg2rad(data_ex.Angle))*repmat(uZvert/100,[1 2])]}; %Z shank coordinates

        if plotIns
            if ~isempty(Fig)
                hold on
            end

            plot3(Coor{2,1},Coor{2,2},Coor{2,3})

        end


    else %%NP1.0

        uX = zeros(2,size(goodU,2));
        uY = zeros(2,size(goodU,2));
        uZ = zeros(2,size(goodU,2));

        if length(cell2mat(data_ex.Matlab3Dcoor))>3

            parts = split(cell2mat(data_ex.Matlab3Dcoor), ',');

            Zstl = surfaceZ(str2double(parts{1}),str2double(parts{2}), brain);

        else

            parts{1} = data_ex.coor_X/100;
            parts{2} = data_ex.coor_Y/100;
            Zstl = surfaceZ(str2double(parts{1}),str2double(parts{2}), brain);

        end

        goodUdepth = NP.chLayoutPositions(2,goodU(1,:));

        verticalDepth = data_ex.Depth -goodUdepth; %depth of unit along vertical axis

        uY(2,:) = (str2double(parts{2})*100) + cos(deg2rad(data_ex.Angle)).*(verticalDepth); %uM
        uY(1,:) = str2double(parts{2}) + cos(deg2rad(data_ex.Angle)).*(verticalDepth/100); %STL

        uZ(2,:) = sin(deg2rad(data_ex.Angle)).*(verticalDepth); %uM
        uZ(1,:) =Zstl-sin(deg2rad(data_ex.Angle)).*(verticalDepth/100); %STL

        uX(2,:) = str2double(parts{1})*100; %uM
        uX(1,:) = str2double(parts{1}); %STL

        Coor = {uX,uY,uZ;... %unit coordinates
            [str2double(parts{1});str2double(parts{1})],... %X shank coordinates STL
            [str2double(parts{2});str2double(parts{2})+cos(deg2rad(data_ex.Angle))*(data_ex.Depth/100)],... %Y shank coordinates
            [Zstl;Zstl - sin(deg2rad(data_ex.Angle))*(data_ex.Depth/100)]}; %Z shank coordinates


        if plotIns
            if ~isempty(Fig)
                hold on
            end

            plot3(Coor{2,1},Coor{2,2},Coor{2,3},'LineWidth',2)

        end
        
        if plotN
            if ~isempty(Fig)
                hold on
            end

            for u = 1:size(goodU,2)

                scatter3(uX(1,u),uY(1,u),uZ(1,u),25,colors(u,:),'filled') 

            end

              view(sideview)
        end




    end

    save(sprintf('Insertion_coordinates-%s',NP.recordingName),'Coor')

else

    Coor = load(sprintf('Insertion_coordinates-%s.mat',NP.recordingName)).Coor;

end

end




%%% Function to get surface point given x and y on 3d model.
function [closest_z]= surfaceZ(x_target,y_target,brain)

vertices = brain.Points;

% Calculate the squared distance between vertices and the target point
distances_squared = (vertices(:, 1) - x_target).^2 + (vertices(:, 2) - y_target).^2;

% Find vertices with distances within a small tolerance
tolerance = 1;  % Adjust the tolerance as needed
close_indices = find(distances_squared < tolerance);

if isempty(close_indices)
    % No vertices are very close
    [~, min_index] = min(distances_squared);
else
    % Among close vertices, find the one with the highest Z coordinate
    [~, max_z_index] = max(vertices(close_indices, 3));
    min_index = close_indices(max_z_index);
end

% Get the Z coordinate of the closest vertex
closest_z = vertices(min_index, 3);


end


