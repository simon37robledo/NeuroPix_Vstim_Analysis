% ----------- Setup ------------
videoFile = 'PV132_3_Camera1_20250505-154817DLC_Resnet50_IR_Eye_movementsApr29shuffle1_snapshot_200_p60_labeled.mp4';
filteredCsv = 'filtered_output.csv';
outputVideo = 'filtered_annotated_video.mp4';

% Read filtered table
T_filtered = readtable(filteredCsv, ...
    'ReadVariableNames', true, ...
    'NumHeaderLines', 3);

% Read DLC header rows
fid = fopen(filteredCsv);
scorerLine = fgetl(fid);
bodypartsLine = fgetl(fid);
coordsLine = fgetl(fid);
fclose(fid);

% Parse body parts
bodyparts = strsplit(bodypartsLine, ',');
coords = strsplit(coordsLine, ',');

% Open video
v = VideoReader(videoFile);
writer = VideoWriter(outputVideo, 'MPEG-4');
writer.FrameRate = v.FrameRate;
open(writer);

% Colormap for drawing (optional)
colors = lines(numel(unique(bodyparts)) / 3);

% ----------- Frame Loop ------------
frameIdx = 1;
while hasFrame(v) && frameIdx <= height(T_filtered)
    frame = readFrame(v);

    % Draw filtered annotations
    imshow(frame); hold on;

    for i = 1:3:length(bodyparts)
        part = strrep(bodyparts{i}, '"', '');  % remove quotes if present
        x_col = i;
        y_col = i + 1;
        likelihood_col = i + 2;

        x = T_filtered{frameIdx, x_col};
        y = T_filtered{frameIdx, y_col};
        p = T_filtered{frameIdx, likelihood_col};

        % Only plot if confidence is high
        if p > 0.8 && all(~isnan([x, y]))
            plot(x, y, 'o', ...
                'MarkerSize', 6, ...
                'MarkerFaceColor', colors(floor(i/3)+1,:), ...
                'MarkerEdgeColor', 'k');
        end
    end

    % Capture frame and write
    F = getframe(gca);
    writeVideo(writer, F.cdata);

    hold off;
    frameIdx = frameIdx + 1;
end

close(writer);
