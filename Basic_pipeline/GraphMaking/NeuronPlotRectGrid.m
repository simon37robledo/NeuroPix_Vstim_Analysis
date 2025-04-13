
coorRect = [0           0        1920        1080];

reduceFactor =10;
coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));
%reduceFactor = min([20 min(sizes)]); %has to be bigger than the smallest ball size
redCoorX = round(coorRect(3)/reduceFactor);
redCoorY = round(coorRect(4)/reduceFactor);

eye_to_monitor_distance = 21.5; % Distance from eye to monitor in cm
pixel_size = 33/(1080/reduceFactor); % Size of one pixel in cm (e.g., 25 micrometers)
monitor_resolution = [redCoorX, redCoorY]; % Width and height in pixels
[theta_x,theta_y] = pixels2eyeDegrees(eye_to_monitor_distance,pixel_size,monitor_resolution);

coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));
redCoorX = round(coorRect(3)/reduceFactor);
redCoorY = round(coorRect(4)/reduceFactor);

theta_x = theta_x(:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY);


RFu = RFrectGrid;

%%%Filter with gaussian:

TwoDGaussian = fspecial('gaussian',floor(size(RFu,2)/(9/2)),size(RFu,1)/9); %increase size of gaussian by 100%.
%RFuSTDirFilt = zeros(size(RFuSTDir));

RFuFilt= zeros(size(RFu));

for ui =1:size(RFu,3) %units

    slice = squeeze(RFu(:,:,ui));

    slicek = conv2(slice,TwoDGaussian,'same');

    RFuFilt(:,:,ui) =slicek;
end


u= 60;
figRF=figure;
imagesc((squeeze(RFuFilt(:,:,u))).*(1000/duration));



c = colorbar;
title(c,'spk/s')

%caxis([0.81 14.59]);

colormap('turbo')
title(sprintf('u-%d',u))

xt = xticks;
xt = xt((1:2:numel(xt)));
xticks(xt);
xticklabels(round(theta_x(1,xt)))

yt = yticks;
yt = yt(1:2:numel(yt));
yticks(yt);
yticklabels(round(theta_y(yt,1)))
%                     xlabel('X degrees')
%                     ylabel('Y degrees')


axis equal tight

saveDir = '\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure';
cd(saveDir)
print(gcf, sprintf('%s-NEM-SameAxisasMB-rectGrid-receptiveField-u-%d.pdf',NP.recordingName,u), '-dpdf', '-r300', '-vector')
