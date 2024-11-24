
data_ex = data(ex,:);

plotN = 1;

plotIns =0;

%function [Coor] = neuronLocation(NP,data_ex,goodU,plotN)

path = convertStringsToChars(string(data_ex.Base_path)+filesep+string(data_ex.Exp_name)+filesep+"Insertion"+string(data_ex.Insertion)...
    +filesep+"catgt_"+string(data_ex.Exp_name)+"_"+string(data_ex.Insertion)+"_g0");
try %%In case it is not run in Vstim computer, which has drives mapped differently
    cd(path)
catch
    originP = cell2mat(extractBetween(path,"\\","\Large_scale"));
    if strcmp(originP,'sil3\data')
        path = replaceBetween(path,"","\Large_scale","W:");
    else
        path = replaceBetween(path,"","\Large_scale","Y:");
    end
    cd(path)
end

if ~exist(string(NP.recordingDir)+filesep+string(NP.recordingName)+"_g0_tcat.imec0.ap_kilosortChanMap.mat")

    SGLXMetaToCoords(string(NP.recordingDir)+filesep+string(NP.recordingName)+"_g0_tcat.imec0.ap.meta");

end

chanMap = load(string(NP.recordingName)+"_g0_tcat.imec0.ap_kilosortChanMap.mat");

if plotN
    path_Brain = '\\132.66.45.127\data\Large_scale_mapping_NP\Model_brain\PV_103_imaged_inclined_45TA.stl';

    brain = stlread(path_Brain);

    figure(1);
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
    view(sideview)

    Ylims = ylim;
    Xlims = xlim;
    Zlims = zlim;

    ymax = 80;
    zmin =-65;

    ylim([Ylims(1)+15 ymax])
    zlim([zmin Zlims(2)])

    % yticks([5:5:ymax-30])
    % zticks([zmin:5:round(Zlims(2))])
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


    %%%Load coordinates from excel data

    if length(cell2mat(data_ex.Matlab3Dcoor))>3 %%Check if the insertion has matlab coordinates

        if string(cell2mat(data_ex.Probe)) == "NP2.0"
            
            parts = split(cell2mat(data_ex.Matlab3Dcoor), ';');

            Xs = split(parts{1},',');
            x1 = str2double(Xs{1});
            x2 = str2double(Xs{2});

            Ys =  split(parts{2},',');
            y1 = str2double(Ys{1});
            y2 = str2double(Ys{2});

            % Vectors for AB and AC
            AB = [0, y2 - y1];
            AC = [x2 - x1, y2 - y1];

            % Calculate the dot product and magnitudes of AB and AC
            dotProduct = dot(AB, AC);
            magAB = norm(AB);
            magAC = norm(AC);

            % Calculate the angle in radians and convert to degrees
            angleRad = acos(dotProduct / (magAB * magAC));

            coorS1 = [x1,y1,surfaceZ(x1,y1, brain)];
            coorS4 = [x1 - 750/100 * cos(angleRad), y1 - 750/100 * sin(angleRad), surfaceZ(x1 - 750/100 * cos(angleRad), y1 - 750/100 * sin(angleRad),brain)];
            coorS2 = [x1 - 750/3/100 * cos(angleRad), y1 - 750/3/100 * sin(angleRad), surfaceZ(x1 - 750/3/100 * cos(angleRad), y1 - 750/3/100 * sin(angleRad),brain)];
            coorS3 = [x1 - 750/(3/2)/100 * cos(angleRad), y1 - 750/(3/2)/100 * sin(angleRad), surfaceZ(x1 - 750/(3/2)/100 * cos(angleRad), y1 - 750/(3/2)/100 * sin(angleRad),brain)];
            
            x = [coorS1(1),coorS2(1),coorS3(1),coorS4(1)];
            y = [coorS1(2),coorS2(2),coorS3(2),coorS4(2)];
            z = [coorS1(3),coorS2(3),coorS3(3),coorS4(3)];

            if plotIns
                scatter3(x,y,z,'bo', 'filled')
            end

            colors = ['bo','ro'];
            j =1;
            for u = 1:length(goodU)


                %Main unit channel position

                %change channel 0 to 1

                ch = goodU(1,34);

                uX = x(chanMap.kcoords(ch));

                uYvert = y(chanMap.kcoords(ch));

                uZvert = max(z) - data_ex.coor_Z/100 - chanMap.ycoords(ch)/100; %I take the max Z because depth is counted from the insertion of first shank

                uZ = sin(deg2rad(data_ex.Angle))*uZvert; %depth of unit along vertical axis

                uY = uYvert-cos(deg2rad(data_ex.Angle))*uZvert;

                scatter3(uX,uY,uZ,colors(j), 'filled')
                view(sideview)

                

                j=j+1;


            end

        else

        end

    end


%     figure;
%     plot()

%      scatter3(animals{a}(:,1),animals{a}(:,2),animals{a}(:,3), 'bo', 'filled', 'MarkerEdgeColor','w', 'MarkerFaceColor',colorA{a}, 'LineWidth',2, 'SizeData',sizeDot);
%         text(animals{a}(:,1)+offsetTextY,animals{a}(:,2)-offsetTextX,animals{a}(:,3)+offsetTextZ, string((1:length(animals{a}))+in),...
%             'Color', 'k','FontSize',sizeText);

end

length(cell2mat(data.Matlab3Dcoor(40)))

%end



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


