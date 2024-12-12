%% Do the sumary graphs. Boxplot for moving ball. 

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

examplesSDG =[1 2 3 4 5 6 7 8 9 10 11 12 13 14 29 30 31 32 40 41 42 43];

%%
animalN = 4;
i=1;
OSI =cell(1,length(examplesSDG));
DSI = cell(1,length(examplesSDG));
PreferAngle = cell(1,length(examplesSDG));
animalID = cell(1,animalN);
DepthUnit = cell(1,length(examplesSDG));

for ex = examplesSDG

     path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
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
    NP = NPAPRecording(path);

    NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;

    tuningCurve = load(sprintf('tuningC-%s',NP.recordingName)).tuningCurve;

    uDir = rad2deg(unique(squeeze(NeuronVals(:,:,5))));
    [preferDir dirInd] = max(tuningCurve,[],2);

    preferDir = uDir(dirInd);

    goodNeurons = load(sprintf('pvalTime-%s',NP.recordingName)).pvalTi;

    goodNeurons = find(goodNeurons<0.05);

    p = NP.convertPhySorting2tIc(NP.recordingDir);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');
    goodUdepth = NP.chLayoutPositions(2,goodU(1,:));
    uDepth = sin(deg2rad(data.Angle(ex))).*goodUdepth(goodNeurons);
    verticalDepth = sin(deg2rad(data.Angle(ex))).*(data.Depth(ex)); %depth of unit along vertical axis
    verticalDepthU= verticalDepth- uDepth;

    OSIi = load(sprintf('Orientation-Tuning-Index-%s',NP.recordingName)).L';
    DSIi = load(sprintf('Direction-Selectivity-Index-%s',NP.recordingName)).DSI';

    %theta{i} = load(sprintf('Angle-prefer-%s',NP.recordingName)).Theta;
    OSI{i} = OSIi(goodNeurons);
    DSI{i} = DSIi(goodNeurons);
    PreferAngle{i} = preferDir(goodNeurons);
    animalID{i} = data.Animal_ID(ex);
    DepthUnit{i} = verticalDepthU;
    i=i+1;

end
%%

group_names = {'OSI', 'DSI'};

h = figure;
daviolinplot({[OSI{:}]' [DSI{:}]'} ,'boxcolors','k','outliers',0,'whiskers',1,...
    'box',2,'boxwidth',1.2,'scatter',1,'scattersize',10,'jitter',1,'violinwidth',1,'jitterspacing',0.1,...
    'xtlabels', {group_names{1} group_names{2}},'linkline',1);
%title(sprintf('Diff-SDG-%s',strrep(NP.recordingName,'_','-')))
ylabel('Orientation Index - Selectivity Index')
set(h,'Color','w');%yline(-10:2.5:15,'LineWidth',0.1,'Alpha',0.3);%ylim([-10,15]);
cd(NP.recordingDir+"\Figs")
print(h, sprintf('SDG-%s.png',NP.recordingName),'-dpng');

%% Swarm plots 

values = [OSI{:} DSI{:}];

[uniqueStrings,~,indexAnimal] = unique([animalID{:}]);
color1 =[];

for i = 1:length(animalID)
    
    color1 = [color1 zeros(1,length(OSI{i}))+indexAnimal(i)];

end

color = [color1 color1];

cats = categorical([ones(1,length(cell2mat(OSI))) ones(1,length(cell2mat(OSI)))+1]);

xticklabels({'OSI','DSI'})


T = table(cats',values',color','VariableNames',{'cats','values','color'});

figure;
colors = lines(3);
swarmchart(T.cats, T.values, 10, T.color, 'filled');


% hold on;
% uniqueGroups = unique(T.color);  % Get unique group values
%  % Get the same colormap used in the plot
% 
% for i = [1 3]
%     % Plot invisible points for legend creation
%     scatter(nan, nan, 100, colors(i,:), 'filled');
%     
% end
% 
% 
% 
hold on
for i =1:2

    groupMean = mean(T.values(T.cats==categorical(i))); 

     plot([i-0.3, i+0.3], [groupMean, groupMean], 'g-', 'LineWidth', 2);  % Short horizontal line for the mean

end
% 
% legend({'PV35', 'PV103', 'SA5'}, 'Location', 'best');
grid on
set(gcf,'Color','w');%
xticklabels({'OSI','DSI'})
ylabel('Score (0-1)')


%% Plot top DSI with lines
rng(42)
DSIc = [DSI{:};ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))];

mean(ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}])))

OSIc = [OSI{:};ones(1,length([DSI{:}]))+1+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))];

indTOP_DSI= find(DSIc(1,:)>prctile(DSIc(1,:),90));

figure;tiledlayout(1,2,"TileSpacing","compact")

nexttile

% Scatter plot for the first group of points
scatter(DSIc(2,indTOP_DSI), DSIc(1,indTOP_DSI), 50, 'b', 'filled');
hold on;
% Scatter plot for the second group of points
scatter(OSIc(2,indTOP_DSI), OSIc(1,indTOP_DSI), 50, 'r', 'filled');
xticks([mean(ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))),mean(ones(1,length([DSI{:}]))+1+0.1+(0.2-0.1)*rand(1,length([DSI{:}])))])
xticklabels([])

plot([DSIc(2,indTOP_DSI);OSIc(2,indTOP_DSI)],[DSIc(1,indTOP_DSI);OSIc(1,indTOP_DSI)],'Color','k');

ylabel('Score (0-1)');
title('Top 90% percentile DSI');
xlim([1 2.3])
grid on;

% Plot top OSI with lines
rng(42)
DSIc = [DSI{:};ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))];

mean(ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}])))

OSIc = [OSI{:};ones(1,length([DSI{:}]))+1+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))];

indTOP_OSI= find(OSIc(1,:)>prctile(OSIc(1,:),90));

nexttile;

% Scatter plot for the first group of points
scatter(DSIc(2,indTOP_OSI), DSIc(1,indTOP_OSI), 50, 'b', 'filled');
hold on;
% Scatter plot for the second group of points
scatter(OSIc(2,indTOP_OSI), OSIc(1,indTOP_OSI), 50, 'r', 'filled');
xticks([mean(ones(1,length([DSI{:}]))+0.1+(0.2-0.1)*rand(1,length([DSI{:}]))),mean(ones(1,length([DSI{:}]))+1+0.1+(0.2-0.1)*rand(1,length([DSI{:}])))])
xticklabels([])

plot([DSIc(2,indTOP_OSI);OSIc(2,indTOP_OSI)],[DSIc(1,indTOP_OSI);OSIc(1,indTOP_OSI)],'Color','k');
title('Top 90% percentile OSI');
legend('DSI','OSI')
xlim([1 2.3])
grid on;
set(gcf,'Color','w')


%% Plot distribution of prefered angles


pAngles = cell2mat(PreferAngle');

DSIv = cell2mat(DSI);
OSIv = cell2mat(OSI);
%pTop = pAngles(DSIv)>0.7);

indTOP_DSI= find(DSIv>prctile(DSIv,75));
pTopDSI = pAngles(indTOP_DSI);

indTOP_OSI= find(OSIv>prctile(OSIv,75));
pTopOSI = pAngles(indTOP_OSI);

figure;
tiledlayout(2,1,"TileSpacing",'compact')
nexttile
binEdges = [-45,45,135,225,315];
histogram(pTopDSI,'binEdges',binEdges)
xticks([0,90,180,270])
ylabel('No. units with DSI > 75pth')
set(gca, 'XTick', []);
nexttile
histogram(pTopOSI,'binEdges',binEdges)
xticks([0,90,180,270])
ylabel('No. units with OSI > 75pth')
set(gcf,'Color','w')

%% Swarm plot for angles
indTOP_DSI75= find(DSIv>prctile(DSIv,75));
indTOP_DSI50= find(DSIv>prctile(DSIv,50));
indTOP_DSI25= find(DSIv>prctile(DSIv,25));

pths = zeros(length(OSIv),1);
pths(indTOP_DSI25) = 3;
pths(indTOP_DSI50) = 2;
pths(indTOP_DSI75) = 1;

T = table(DSIv',OSIv',pAngles,color1',pths, cell2mat(DepthUnit)','VariableNames',{'DSI','OSI','pAngles','Animals','percentiles','depth'});

%1 animal
% figure;
% swarmchart(T.pAngles, T.DSI, 10, T.percentiles, 'filled');
% xticks([0:45:315])
% xlabel('angle')
% ylabel('DSI score')

%%  Manual violin plot per angle
uniqueAngles = 0:45:315;

groups = T.DSI;
group_names = {'0','45','90','135','180','215','270','315'};
for a = 1:8
    groupA{a} = groups(T.pAngles == uniqueAngles(a));
end


% Sample data (1x8 cell array with 8 groups of data)
data = groupA;

% Create figure
figure;
hold on;

% Define colors for the 8 groups (optional: customize this)
colors = lines(8);  % Use a colormap with 8 distinct colors

% Find the global min and max across all groups
globalMin = min(cellfun(@min, data));  % Global minimum across all data
globalMax = max(cellfun(@max, data));  % Global maximum across all data

% Number of bins for the histograms
nBins = 10;

binEdges = linspace(globalMin-0.05, globalMax+0.05, nBins + 1);  % Consistent bin edges for all groups

% Define a spacing factor to separate the violin plots more (adjust this as needed)
spacingFactor = 6;  % Increase this to create more separation between groups

% Loop over each group and create a smoothed, mirrored histogram
lineWidth = 2;
for i = 1:length(data)
    % Get the current group's data
    currentData = data{i};
    
    % Compute histogram counts using consistent bin edges
    [counts, ~] = histcounts(currentData, binEdges, 'Normalization', 'pdf'); 
       % Get bin centers from the consistent bin edges
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;    
    % Interpolate to smooth the counts
    binInterp = linspace(binCenters(1), binCenters(end), 100);  % More points for smoothing
    countsInterp = interp1(binCenters, counts, binInterp, 'spline');  % Smooth using spline interpolation
    
    % Plot the violin shape using fill by mirroring the density and adding more spacing
    fill([-countsInterp, fliplr(countsInterp)] + i * spacingFactor, ...
        [binInterp, fliplr(binInterp)], colors(i, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
    % Calculate the 25th, 50th (median), and 75th percentiles
    lowerBound = prctile(currentData, 25);
    upperBound = prctile(currentData, 75);
    medianValue = median(currentData);
    meanValue = mean(currentData);
    
    % Plot the boxplot (25% to 75% range)
    line([i*spacingFactor - lineWidth, i*spacingFactor + lineWidth], [lowerBound, lowerBound], 'Color', 'k', 'LineWidth', 2);  % Lower bound
    line([i*spacingFactor - lineWidth, i*spacingFactor + lineWidth], [upperBound, upperBound], 'Color', 'k', 'LineWidth', 2);  % Upper bound
    line([i*spacingFactor, i*spacingFactor], [lowerBound, upperBound], 'Color', 'k', 'LineWidth', 2);  % Vertical line connecting bounds
    
    % Plot the median line
    %line([i*spacingFactor - 1, i*spacingFactor + 1], [medianValue, medianValue], 'Color', 'k', 'LineWidth', 2);  % Median line
    % Plot the mean line
    line([i*spacingFactor - 1, i*spacingFactor + 1], [meanValue, meanValue], 'Color', 'k', 'LineWidth', 2);  % Median line
    
    % Scatter all data points, randomly within the group's x-axis position
    jitterX = (rand(size(currentData)) - 0.5) * 2;  % Random jitter for x-axis, scaled by 0.3
    scatter(i * spacingFactor + jitterX, currentData, 10, 'k', 'filled','MarkerFaceAlpha',0.3);  % Scatter all data points
end

% Customize the plot
xlim([0.5, length(data) * spacingFactor + 3]);  % Set x-axis limits to show all groups
xticks(spacingFactor * (1:length(data)));  % Set x-axis ticks at each spaced group position
xticklabels(group_names);  % Label each group
ylabel('DSI score');
xlabel('Angles');
%title('Smoothed Histogram-Based Violin Plot with Boxplot and Scatter');
grid on;
hold off;
set(gcf,'Color','w')

%% Depth plots per angle

%%  Manual violin plot per angle
uniqueAngles = 0:45:315;

groups = T.depth;
group_names = {'0','45','90','135','180','215','270','315'};
pthcond = T.DSI > prctile(DSIv,50);
for a = 1:8
    groupA{a} = groups(T.pAngles == uniqueAngles(a) & groups>=0 & T.Animals==4);
end


% Sample data (1x8 cell array with 8 groups of data)
data = groupA;

% Create figure
figure;
hold on;

% Define colors for the 8 groups (optional: customize this)
colors = lines(8);  % Use a colormap with 8 distinct colors

% Find the global min and max across all groups
globalMin = min(cellfun(@min, data));  % Global minimum across all data
globalMax = max(cellfun(@max, data));  % Global maximum across all data

% Number of bins for the histograms
nBins = 10;

binEdges = linspace(0-250, 4000+250, nBins + 1);  % Consistent bin edges for all groups

% Define a spacing factor to separate the violin plots more (adjust this as needed)
spacingFactor = 0.002;  % Increase this to create more separation between groups

% Loop over each group and create a smoothed, mirrored histogram
lineWidth = 0.0005;
for i = 1:length(data)
    % Get the current group's data
    currentData = data{i};
    
    % Compute histogram counts using consistent bin edges
    [counts, ~] = histcounts(currentData, binEdges, 'Normalization', 'pdf'); 
       % Get bin centers from the consistent bin edges
    binCenters = (binEdges(1:end-1) + binEdges(2:end)) / 2;    
    % Interpolate to smooth the counts
    binInterp = linspace(binCenters(1), binCenters(end), 100);  % More points for smoothing
    countsInterp = interp1(binCenters, counts, binInterp, 'spline');  % Smooth using spline interpolation
    
    % Plot the violin shape using fill by mirroring the density and adding more spacing
    fill([-countsInterp, fliplr(countsInterp)] + i * spacingFactor, ...
        [binInterp, fliplr(binInterp)], colors(i, :), 'FaceAlpha', 0.5, 'EdgeColor', 'none');
    
    % Calculate the 25th, 50th (median), and 75th percentiles
    lowerBound = prctile(currentData, 25);
    upperBound = prctile(currentData, 75);
    medianValue = median(currentData);
    meanValue = mean(currentData);
    
    % Plot the boxplot (25% to 75% range)
    line([i*spacingFactor - lineWidth, i*spacingFactor + lineWidth], [lowerBound, lowerBound], 'Color', 'k', 'LineWidth', 2);  % Lower bound
    line([i*spacingFactor - lineWidth, i*spacingFactor + lineWidth], [upperBound, upperBound], 'Color', 'k', 'LineWidth', 2);  % Upper bound
    line([i*spacingFactor, i*spacingFactor], [lowerBound, upperBound], 'Color', 'k', 'LineWidth', 2);  % Vertical line connecting bounds
    
    % Plot the median line
    %line([i*spacingFactor - 1, i*spacingFactor + 1], [medianValue, medianValue], 'Color', 'k', 'LineWidth', 2);  % Median line
    % Plot the mean line
    line([i*spacingFactor - lineWidth, i*spacingFactor + lineWidth], [meanValue, meanValue], 'Color', 'k', 'LineWidth', 2);  % Median line
    
    % Scatter all data points, randomly within the group's x-axis position
    jitterX = (rand(size(currentData)) - 0.5) * lineWidth;  % Random jitter for x-axis, scaled by 0.3
    scatter(i * spacingFactor + jitterX, currentData, 10, 'k', 'filled','MarkerFaceAlpha',0.3);  % Scatter all data points
end

% Customize the plot
xlim([-spacingFactor, length(data) * spacingFactor]+spacingFactor);  % Set x-axis limits to show all groups
ylim([0 4000])
xticks(spacingFactor * (1:length(data)));  % Set x-axis ticks at each spaced group position
xticklabels(group_names);  % Label each group
ylabel('Depth (um)');
xlabel('Angles');
%title('Smoothed Histogram-Based Violin Plot with Boxplot and Scatter');
grid on;
hold off;
set(gcf,'Color','w')






















