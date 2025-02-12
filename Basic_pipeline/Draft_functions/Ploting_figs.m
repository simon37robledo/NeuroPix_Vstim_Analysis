%Plotting figures of saved data:
RfuMov = load('RFu_MovingBall-PV139_Experiment_6_2_24_1.mat').normRFu;

figure;imagesc(squeeze(RfuMov(1,:,42:end-42,51)));colorbar;clim([caxisLims(:,j)]);figure;imagesc(normRFu(:,:,51));colorbar
set(gcf, 'Color', 'w');
xlabel('X screen (10 pixels)');ylabel('Y screen (10 pixels)')

%%
RfuDir = load('RFuDirec_MovingBall-PV139_Experiment_6_2_24_1.mat').DirecSum;
for i =1:size(RfuDir,1)
figure;imagesc(squeeze(RfuDir(i,:,42:end-42,65)));cb = colorbar; cb.Label.String = 'Response/Baseline';title(sprintf('Theta = %d',90*(i-1)));set(gcf, 'Color', 'w');
xlabel('X screen (10 pixels)');ylabel('Y screen (10 pixels)');%caxis([-1 2.5]);
end

RFuConv = load('RFu-51.mat').RFu;

implay(RFuConv)

NormM = load('NormMatrix_MovingBall-PV139_Experiment_6_2_24_1.mat').normMatrix;

figure;imagesc(squeeze(mean(NormM(:,:,42:end-42,65))));cb = colorbar; cb.Label.String = 'Occupancy*Baseline';title(sprintf('r = %d',5*i));set(gcf, 'Color', 'w');
xlabel('X screen (10 pixels)');ylabel('Y screen (10 pixels)');

%% Comparison RG, MB

selectU = [51,65,21,26,32,34,42,47,53,58];

caxisLims = zeros(2,length(selectU));

j =1;
figure;
t= tiledlayout(2,length(selectU),"TileSpacing","compact");
for i =selectU%size(RfuDir,1
    nexttile
    imagesc(squeeze(mean(RfuMov(:,:,42:end-42,i))));title(sprintf('Unit %d',i));set(gcf, 'Color', 'w');
    caxisLims(:,j) = [min(squeeze(mean(RfuMov(:,:,42:end-42,i))),[],'all') max(squeeze(mean(RfuMov(:,:,42:end-42,i))),[],'all')];
    j = j+1;
end

j=1;
for i =selectU%size(RfuDir,1)
    nexttile;
    imagesc(flip(flip(squeeze(normRFu(:,:,i)),2),1));%cb = colorbar;set(gcf, 'Color', 'w');
    clim([caxisLims(:,j)]);
    j = j+1;
end
% Get all axes in the tiled layout
axs = findall(t,'Type','axes');

% Loop through each axes and remove y ticks
for i = 1:length(axs)
    axs(i).YTick = [];
    axs(i).XTick = [];
end

xlabel('X screen (10 pixels)');ylabel('Y screen (10 pixels)');
cb.Label.String = 'Response/Baseline';

figure;imagesc(squeeze(normRFu(:,:,21)));cb = colorbar; cb.Label.String = 'Response/Baseline';set(gcf, 'Color', 'w');
xlabel('X screen (10 pixels)');ylabel('Y screen (10 pixels)');

%%
u =65;


 [nT,nB2] = size(Mr2);

 fig = figure;

 imagesc(squeeze(Mr2(:,u,:)));colormap(flipud(gray(64)));
 %Plot stim start:
 xline(preBase/bin,'k', LineWidth=1.5)
 %Plot stim end:
 xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5)
 ylabel('Trials');xlabel('Time (ms)');
 title(sprintf('U.%d-R.%.3f-B.%.3f',u,max_mean_value(u),spkRateBM(u)));

 %xticks([0.5 (preBase/bin):10:nB])
 %xticklabels([-preBase 0:10*bin:nB*bin])
 yticklabels([yticks]*mergeTrials)
 %Directions
 v = nT/direcN:nT/direcN:nT-1;
 yline(v+0.5,'r', LineWidth=3);
 %Offsets
 v = nT/(direcN*offsetN):nT/(direcN*offsetN):nT-1;
 yline(v+0.5,'b', LineWidth=2);
 %sizes
 v = nT/(direcN*offsetN*sizeN):nT/(direcN*offsetN*sizeN):nT-1;
 yline(v+0.5, LineWidth=0.5);

 %                             hcb = colorbar();
 %                             title(hcb,'Spikes/sec');
 %caxis([0 max(0.2,max(max_mean_value(u)))])
 hold on
 %Plot rectangle:

 [MaxWin mI] = max(squeeze(NeuronVals(u,:,1)));

 mxTrials = squeeze(NeuronVals(u,mI,2));

 mxBin = squeeze(NeuronVals(u,mI,3));


 rectangle('Position',[mxBin,mxTrials,window_size(2),1],...
     'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
                hold off
                prettify_plot


                %%
tic
duration = 400;
% Convert image to double precision for computation
for u =1:length(goodU)
img = squeeze(Mr2(:,u,:));

% Define window size
windowSize = [trialDiv round(duration/bin)]; % for a 5x5 window

% Define a function to calculate mean
meanFun = @(x) mean(x(:));

% Apply the function to the image using a sliding window
meanImg = nlfilter(img, windowSize, meanFun);

% Find the location of the maximum mean 
[maxMean, maxIndex] = max(meanImg(:));
% Convert the linear index to row, column subscripts
[row, col] = ind2sub(size(meanImg), maxIndex);
end
toc

% Now you have the top-left corner of the window with the highest mean
% You can extract this window from the original image
highestMeanWindow = img(row:row+windowSize(1)-1, col:col+windowSize(2)-1);

%figure;imagesc(meanImg)
