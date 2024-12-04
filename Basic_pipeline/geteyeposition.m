
%% get eye positions for a choosen recording (combine with head movement), gives eye x position for sorted onStim, and in ellipse x-y positions and diameter

%function [meanHead_pos, ellipse] = geteyeposition(animal_Nr, plotting)
% 
% folder = '/media/sil1/Data/Turtle neurotassel recordings/Ephys Video for DLC/eye_labeling_turtles-Milan-2023-10-01/videos';
% folders_all = dir(folder);
% folders_all = {folders_all.name};
% folders_all = folders_all(contains(folders_all, '.csv')); %this will give you the name of all the files in the folder
% animals_alll = squeeze(split(folders_all, {'trial', 'DLC'}));
% 
% for i = 1:length(animals_alll)
%     animals_all(i) = append(animals_alll(i,1), animals_alll(i,2));
%     if contains(animals_all{i}, 'LE')
%         animals_all{i}(end-2:end) = [];
%     end
% end
% if ~exist('plotting', 'var')
%     plotting = 1;
% end

%% start loop -- no loop

% pp = find(contains(animals_all, animal_Nr), 1, 'last'); % finds the recording for the animal --- only LE!!
% 
% %% get turtle folder
% 
% turtle_Nr = split(animal_Nr,'_');
% turtle_Nr = turtle_Nr{2};
% turtle_Nr(1) = 'T';
% turtle_Nr = append(turtle_Nr(1:end-2), ' ', turtle_Nr(end-1:end));
% turtle_folder = fullfile('/media/sil1/Data/Turtle neurotassel recordings', turtle_Nr);
% 
% folder_for_saving = fullfile(turtle_folder, animal_Nr);
% filename = fullfile(folder, string(folders_all(pp))); %gives the full path to the file you want to look at
% 
% if contains(filename, 'LE')
%     eye = 'LE';
% else
%     eye = 'RE';
% end

%% get video path
% videopath = dir(fullfile(folder_for_saving, 'eye_videos', eye));
% videopath1 = {videopath.name};
% videopath1 = videopath1{3};
% videopatt = dir(fullfile(folder_for_saving, 'eye_videos', eye, videopath1));
% videopath2 = {videopatt.name};
% videopath2 = videopath2{contains({videopatt.name}, '.mp4')};
% 
% videopath = fullfile(folder_for_saving, 'eye_videos', eye, videopath1, videopath2);
% %% load things
% 
% behavtable = readtable(filename , 'HeaderLines', 1); %reads the .csv table into matlab
% 
% data_folder = fullfile(folder_for_saving, 'Record Node 101/');
% if ~exist(data_folder, "dir")
%     data_folder = fullfile(folder_for_saving, 'Record Node 113/'); % tries to find the right record node
%     if ~exist(data_folder, "dir")
%         data_folder = fullfile(folder_for_saving, 'Record Node 105/');
%     end
% end
% 
% % dataObj
% dataObj=OERecording(data_folder);
% digiData = dataObj.getTrigger; % gives start end end points of each trigger
% digiData = correct_digiData_offset(digiData, folder_for_saving); %corrects the digiData if there is an offset saved
% %     load(fullfile(folder_for_saving, 'good_channels.mat'), 'good_channels');
% 
% % VS file
% if ~isfolder(fullfile(folder_for_saving,"data/"))
%     %         continue
% end
% stim_path=dir(fullfile(folder_for_saving,"data/")); %lists files in folder
% stim_path1={stim_path.name}; %converts it to a cell array to work with
% vStimFile=fullfile(folder_for_saving,"data/"+string(stim_path1(find(contains(stim_path1, '.csv'),1)))); %takes the first file with a csv extension
% VS=readtable(vStimFile); % imports the excel file from psychopy
% [~,~, stim_type] = accXml(VS);
% 
% % run only on novelty
% if ~contains(stim_type, 'novelty')
%     %         continue
% end
% 
% % onStim
% load(fullfile(dataObj.recordingDir,'analysis','frameShifts'), 'frameShifts'); %frameshifts contains on and off of the photodiode each in one cell
% onStim = frameShifts{1}([1;1+find(diff(frameShifts{1})>2500)]);
% 
% % sort the onStim
% [~, ~, ~, onStim] = onStimPerStimulus(stim_type, onStim, VS); % gives a sorted onStim
% 
% % clean behavtable
% [behavtable] = clean_dlc_table(behavtable);
% 
% %% define points
% % VR = VideoReader(videopath);
% xpoints = [behavtable.Pupil_12, behavtable.Pupil_6, behavtable.Pupil_9, behavtable.Pupil_3, behavtable.Pupil_10, behavtable.Pupil_1, behavtable.Pupil_4, behavtable.Pupil_8];
% ypoints = [behavtable.Pupil_12_1, behavtable.Pupil_6_1, behavtable.Pupil_9_1, behavtable.Pupil_3_1, behavtable.Pupil_10_1, behavtable.Pupil_1_1, behavtable.Pupil_4_1, behavtable.Pupil_8_1];
% 
% 
% %% fit the ellipse for each frame
% warning('off', 'MATLAB:illConditionedMatrix'); %disables warning
% for i=1:length(behavtable.Pupil_12) % 6000%:60000
%     % frame = read(VR,i);
%     % imshow(frame);
%     % hold on
%     % xpoints = [behavtable.Pupil_12(i), behavtable.Pupil_6(i), behavtable.Pupil_9(i), behavtable.Pupil_3(i), behavtable.Pupil_10(i), behavtable.Pupil_1(i),behavtable.Pupil_4(i), behavtable.Pupil_8(i)];
%     % ypoints = [behavtable.Pupil_12_1(i), behavtable.Pupil_6_1(i), behavtable.Pupil_9_1(i), behavtable.Pupil_3_1(i), behavtable.Pupil_10_1(i), behavtable.Pupil_1_1(i),behavtable.Pupil_4_1(i), behavtable.Pupil_8_1(i)];
%     ellips = fit_ellipse(xpoints(i,:), ypoints(i,:));
%     if isempty(ellips) %if there is no ellipse possible
%         ellips.X0_in = NaN;
%         ellips.Y0_in = NaN;
%         ellips.phi = NaN;
%         ellips.a = NaN;
%         ellips.b = NaN;
%     end
%     if isempty(ellips.X0_in)
%         ellips.X0_in = NaN;
%         ellips.Y0_in = NaN;
%         ellips.phi = NaN;
%         ellips.a = NaN;
%         ellips.b = NaN;
%     end
%     ellipse.mid(i,:) = [ellips.X0_in, ellips.Y0_in];
%     ellipse.orient(i,:) = ellips.phi*(180/pi); %changes randomly, use pupil_12 and pupil_6 to calculate angle
%     ellipse.radius(i,:) =  mean([ellips.a, ellips.b]); %sqrt(ellips.a*ellipse.b); %[ellips.a, ellips.b]; % this takes the two radiuses of the ellipse and takes the middle
%     % ellipse.ratio(i,:) = ellips.short_axis/ellips.long_axis;
% 
% 
% 
%     if ~exist('ellipsTali', 'var') %get a variable for talis analysis
%         ellipsTali = ellips;
%     else
%         fns = fieldnames(ellips);
%         for i = 1:length(fns)
%             ellipsTali.(fns{i}) = horzcat(ellipsTali.(fns{i}),ellips.(fns{i}));
%         end
%     end
% 
% end
% 
% ellipse.uncorr_mid = ellipse.mid;
% ellipse.speed = hypot(diff(ellipse.mid(:,1)), diff(ellipse.mid(:,2))); % speed in pixel per frame
% 
% % pause(0.5)
% % for i= 1:60
% %     frame = read(VR,150954+i);
% %         imshow(frame);
% % hold on
% % plot(ellipse.mid(150954+i,1),ellipse.mid(69954+i,2), '*')
% % % end


%% gets rostral and caudal points and their midpoint+diameter (--- maybe adjust this using binning if camera moves during recording)
dots_caud = [median(behavtable.Caudal_edge, 'omitmissing'), median(behavtable.Caudal_edge_1, 'omitmissing')];
dots_rost = [median(behavtable.Rostral_edge, 'omitmissing'), median(behavtable.Rostral_edge_1, 'omitmissing')];


d = pdist([dots_caud; dots_rost],'euclidean');
midpoint_eyesocket = dots_caud+0.5*(dots_rost-dots_caud);
%  plot(midpoint_eyesocket(1), midpoint_eyesocket(2), 'r*')
%  viscircles(midpoint_eyesocket,d/2)
ellipse.uncorr_rost = dots_rost;
ellipse.uncorr_caud = dots_caud;
ellipse.uncorr_mid_socket = midpoint_eyesocket;


%% put here for if manual labeling exists
if exist(append(folder, '/eye_socket_manual_', animal_Nr, '.mat'), 'file')
    load(append(folder, '/eye_socket_manual_', animal_Nr, '.mat'))
    ellips_socket = fit_ellipse(socket_all(:,1), socket_all(:,2));
    eyesocket_middle = [ellips_socket.X0_in, ellips_socket.Y0_in];

    midpoint_eyesocket = eyesocket_middle; % change middle and dots caud and rost
    dots_caud = socket_all(3,:);
    dots_rost = socket_all(4,:);
    d = pdist([dots_caud; dots_rost],'euclidean');
    ellipse.uncorr_mid_socket = midpoint_eyesocket;
else
    warning('No manual labels')
end

%% make points around the middle and roatate so that caudal edge is always the same direction
new_middle = ellipse.mid-midpoint_eyesocket;
new_caud = dots_caud-midpoint_eyesocket;
new_rost = dots_rost-midpoint_eyesocket;


%% calculate eye angle with mean eyemid = eye middle
% first calculate in not roated picture the x and y angle of the mean eye
% mid and of all the eye mid points. 
% then subtract the mean eyemid angles to make it 0,0
% then use the rotation matrix to rotate the x,y angle of each eyemid
% get the new x-angle -- in angles from the meaneyemid



%% x-y angle in not rotated picture
d2 = d*0.85; % because d is a little overestimation
radius_socket = d2/2;

for i = 1:length(new_middle)
    dist(i) = norm(new_middle(i,:)); %distance of each point from middle
end
new_middle(dist>radius_socket, :) = NaN;

% mean eyemid
x_eyemid_mean = mean(new_middle(:, 1), 'omitmissing');
y_eyemid_mean = mean(new_middle(:, 2), 'omitmissing');
% all eyemids
x = new_middle(:, 1);
y = new_middle(:, 2);

% remove outliers that are outside the boundaries - false labels or camera move - should be removed already
% x(x>=radius_socket) = radius_socket;
% x(x<=-radius_socket) = -radius_socket;
% y(y>=radius_socket) = radius_socket;
% y(y<=-radius_socket) = -radius_socket;



% x angle mean
radius = sqrt(radius_socket^2-y_eyemid_mean^2);
theta2 = acos(x_eyemid_mean/radius);
x_mean_angle = (rad2deg(theta2)-90)*-1; % negative angles = eye goes to the Caudal side, positive angles - goes to the rostral
% y angle mean
radius = sqrt(radius_socket^2-x_eyemid_mean^2);
theta2 = acos(y_eyemid_mean/radius);
y_mean_angle = (rad2deg(theta2)-90)*-1; % negative angles = ??

% x angle, pythagoras x^2+ydistance^2 = radius^2
radius = sqrt(radius_socket^2-y.^2);
theta2 = acos(x./radius);
xangle = (rad2deg(theta2)-90)*-1; % negative angles = eye goes to the Caudal side, positive angles - goes to the rostral
% y angle
radius = sqrt(radius_socket^2-x.^2);
theta2 = acos(y./radius);
yangle = (rad2deg(theta2)-90)*-1; % negative angles = ?? up or down?

% subtract mean angle
xangle = xangle-x_mean_angle;
yangle = yangle-y_mean_angle;
xyangle = [xangle, yangle]; % make x y angle a point
%% rotate angle x-y to get resulting x angle

% angle of streak
new9_1 = behavtable.Pupil_9-behavtable.Pupil_3; % make point 9 the 0,0 and calculate angle to point 3 from there
new9_2 = behavtable.Pupil_9_1-behavtable.Pupil_3_1;
theta = 180-atan2d(new9_2, new9_1); % angle to rotate in degrees
ellipse.streakangle = theta;
% theta = zeros(length(theta),1); %removes rotation
% ---------------------------------------rotate in right direction ----------------------------------- works for LE - check RE
for i=1:length(theta)
    R = [cosd(theta(i)) -sind(theta(i)); sind(theta(i)) cosd(theta(i))]; %rotation matrix
    % Rotate the points
    point = xyangle'; % shifted midpoints of the ellipse
    rotangle(:,i) = R*point(:,i);
end

%% delete outliers bigger than 90° -- shouldnt exist anymore

for i = 1:length(rotangle)
    dist(i) = norm(rotangle(:,i)); %distance of each point from middle
end
rotangle(dist>90,:) = NaN;

xangle_rot = rotangle(1,:)';
yangle_rot = rotangle(2,:)';

ellipse.midangle = rotangle';
%% get eye position for each stimulus

nr_frames = 24; % 24 = 0.4s window after stimulus onset

if contains(eye, 'RE')
    digiDataVid = digiData{9}; % for RE
end
if contains(eye, 'LE')
    digiDataVid = digiData{15}; % for LE
end
% check if triggers fit or maybe they are switched
if length(digiDataVid)-length(ellipse.midangle)~=0

    digleng1 = length(digiData{9})-length(ellipse.midangle);
    digleng2 = length(digiData{15})-length(ellipse.midangle);
    warning(append('Digital triggers dont fit by ', string(length(digiDataVid)-length(ellipse.midangle)), 'other trigger fits better?, current eye:', eye, ' - LE:', string(digleng2), ' RE:', string(digleng1)))
end

% get the mean eye position for each stimulus trial
for i=1:length(onStim)
    trigger_frame(i) = find(digiDataVid>(onStim(i)),1)-1; %finds the triggers numbers (=frame number) after stimOn subtract 1 to be before stimon
    meanHead_pos(i, 1) = mean(ellipse.midangle(trigger_frame(i):(trigger_frame(i)+nr_frames),1));
    meanHead_pos(i, 2) = mean(ellipse.midangle(trigger_frame(i):(trigger_frame(i)+nr_frames),2));
    meanSpeed(i,1) = mean(ellipse.speed(trigger_frame(i):(trigger_frame(i)+nr_frames))); %speed after stimon
    meanSpeed(i,2) = mean(ellipse.speed(trigger_frame(i)-(1+nr_frames):(trigger_frame(i)-1))); % speed before stimon
    meanHead_pos_firstframe(i,:) = ellipse.midangle(trigger_frame(i),:);

    % tali_mean(i) = mean(tali_behav.phi(trigger_frame(i):(trigger_frame(i)+nr_frames)));
end

ellipse.meanspeed = meanSpeed; % in °/frame now


%% estimates angle of eye
ellipse.angle = meanHead_pos(:,1);
ellipse.angley = meanHead_pos(:,2);

%% estimate all angles for correlation with head angle
ellipse.angle_all = ellipse.midangle(:,1);

%% firstframe angle
ellipse.angle_firstframe = meanHead_pos_firstframe;

%% figure with eye positions in eyesocket circle
% plotting = 1;
if plotting == 1
    figure
    set(gcf, 'WindowState', 'maximized')
    % plot(meanHead_pos(:,1), meanHead_pos(:,2), '*')
        plot(ellipse.mid(:,1), ellipse.mid(:,2), '*')

    % text(meanHead_pos(:,1), meanHead_pos(:,2),string(round(tali_mean)))
        % text(meanHead_pos(:,1), meanHead_pos(:,2),string(round(ellipse.angle)))

    hold on
    % plot(ellipse.rostral(1), ellipse.rostral(2), 'r*')
    % plot(ellipse.caudal(1), ellipse.caudal(2), 'r*')
    plot(0,0, 'r*')
    viscircles([0,0],0.5) % plot the circle
    xticks([-0.5,0.5])
    xticklabels({'Backward', 'Forward'})
    yticks([-0.5,0.5])
    yticklabels({'Bottom', 'Top'})
    title(append(animals_all{pp}, '_', stim_type), 'Interpreter','none')
    daspect([1 1 1]) % makes the plot a square

    if ~isfolder(fullfile(turtle_folder, 'eyepos'))
        mkdir(fullfile(turtle_folder, 'eyepos'))
    end
    pause(0.2)
    saveas(gca,string(fullfile(turtle_folder, 'eyepos', append('eyepos_new_', animals_all{pp}, '_', stim_type, '.png'))))

end




%% figure to check if eyemid and positions make sense
% x=0;
% figure
% for i=trigger_frame
%     x=x+1;
%     frame = read(VR,i);
%     imshow(frame);
%     hold on
%     plot(meanHead_pos(x), meanHead_pos_1(x), 'r*')
%     plot(dots_rost(1), dots_rost(2), 'c*')
%     plot(dots_caud(1), dots_caud(2), 'c*')
%     viscircles([midpoint_eyesocket(1), midpoint_eyesocket(2)],d/2)
%     plot(midpoint_eyesocket(1), midpoint_eyesocket(2), 'c*')
%     pause(0.1)
%end