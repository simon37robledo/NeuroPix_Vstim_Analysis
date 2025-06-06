%%%RF analysis


cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

indexSummary =[];

locationRF = [];

in =1;
totalN = 0;

perDirection =0;

for ex = [51] 
    %%%%%%%%%%%% Load data and data paremeters

    %1. Load NP class
    path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    try %%In case it is not run in Vstim computer, which has drives mapped differently
        cd(pathE)
    catch
        try
            originP = cell2mat(extractBetween(path,"\\","\Large_scale"));
            if strcmp(originP,'sil3\data')
                path = replaceBetween(path,"","\Large_scale","W:");
            else
                path = replaceBetween(path,"","\Large_scale","Y:");
            end

            cd(path)
        catch

            if strcmp(originP,'sil3\data')
                path = replaceBetween(path,"","\Large_scale","\\sil3\data");
            else
                path = replaceBetween(path,"","\Large_scale","\\sil1\data");
            end
            cd(path)

        end
    end
    NP = NPAPRecording(path);

    p = NP.convertPhySorting2tIc(NP.recordingDir);

    label = string(p.label');
    goodU = p.ic(:,label == 'good');

    totalN = totalN+size(goodU,2);



    N_bootstrap =1000;

    cd(NP.recordingDir)
    respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;

   
    sign = 0.05; %%%Significance level used to calculate receptive fields
    respU = find(respNeuronsMB<sign);

% %     NeuronPlotMovingBall(data,ex,respU,...
% %     'savePlot',1,'saveDir','W:\Large_scale_mapping_NP\Figs paper\receptiveFieldExamples','ReField',0,'noEyeMoves',1,'DivisionType', 'XY','EyeMovements',1)
% 
%     %%noEyeMoves
% 
%     NMRF = squeeze(load(sprintf('NEM-RFuSTDirSizeFilt-Q1-Div-X-%s',NP.recordingName)).RFuSTDirSizeFilt);
% 
% 
%     %%EyeMoves
% 
%     RF = squeeze(load(sprintf('RFuSTDirSizeFilt-%s',NP.recordingName)).RFuSTDirSizeFilt);
% 
%     nX = size(RF,3);
%     nY = size(RF,2);
% 
%     index = zeros(4,numel(respU)*4);
%     %Rows: 1. indexvalue NEM; 2. value
% 
%     indexFinal = zeros(5,numel(respU)*4);
%     ui =1;
% 
%     for u = 1:numel(respU) %%signifficance should be 0.05 (default for receptive field calculation
%          
% 
% 
%             NMRFi = abs(squeeze(NMRF(:,:,1+(nX-nY)/2:nY+(nX-nY)/2,u)));
%             RFi = abs(squeeze(RF(:,:,1+(nX-nY)/2:nY+(nX-nY)/2,u)));
% 
%             maxAllNEM = -inf;
%             maxAllEM = -inf;
% 
%             for d = 1:4 %directions
% 
%                 %%%%NEM%%%%%%%%%%%%%%%%%%%
%       
%                 A = squeeze(NMRFi(d,:,:));
%                 % Define the kernel size
%                 kernelSize = 12;
% 
%                 % Compute the 2D moving mean using convolution
%                 meanA = conv2(A, ones(kernelSize) / (kernelSize^2), 'valid');
% 
%                 % Find the maximum mean value and its position
%                 [maxMean, linearIdx] = max(meanA(:));
%                 [rowMax, colMax] = ind2sub(size(meanA), linearIdx);
% 
%                 % Convert to original matrix indices
%                 rowStart = rowMax;
%                 rowEnd = rowMax + kernelSize - 1;
%                 colStart = colMax;
%                 colEnd = colMax + kernelSize - 1;
% 
%                 % Create a mask for the maximal kernel
%                 mask = true(size(A));  % Start with all elements as true
%                 mask(rowStart:rowEnd, colStart:colEnd) = false; % Exclude max kernel region
% 
%                 % Compute the mean of all values outside the max kernel
% 
%                 Aout = A(mask);
% 
%                 %[maxMean] = max(Aout(:));
% 
%                 outsideMean = mean(A(mask));
% 
%                 %%%Calculate centersX of max receptive field
% 
%                 
% 
%                 %%%Calcutate centersY
% 
%                 %%%Calculate locality index           
% 
%                 %index(i,ui) = maxMean/max(meanA,[],'all') - outsideMean/max(meanA,[],'all');
% 
%                 index(1,ui) = maxMean;
%                 index(2,ui) = outsideMean;
% 
%                 if maxMean>maxAllNEM
%                     maxAllNEM = maxMean;
%                 end
% 
%                 %%%%%%%%%%%%% EM %%%%%%%%%%%%
% 
%                 A = squeeze(RFi(d,:,:));
%                 % Define the kernel size
%                 kernelSize = 12;
% 
%                 % Compute the 2D moving mean using convolution
%                 meanA = conv2(A, ones(kernelSize) / (kernelSize^2), 'valid');
% 
%                 % Find the maximum mean value and its position
%                 [maxMean, linearIdx] = max(meanA(:));
%                 [rowMax, colMax] = ind2sub(size(meanA), linearIdx);
% 
%                 % Convert to original matrix indices
%                 rowStart = rowMax;
%                 rowEnd = rowMax + kernelSize - 1;
%                 colStart = colMax;
%                 colEnd = colMax + kernelSize - 1;
% 
%                 % Create a mask for the maximal kernel
%                 mask = true(size(A));  % Start with all elements as true
%                 mask(rowStart:rowEnd, colStart:colEnd) = false; % Exclude max kernel region
% 
%                 % Compute the mean of all values outside the max kernel
% 
%                 Aout = A(mask);
% 
%                 %[maxMean] = max(Aout(:));
% 
%                 outsideMean = mean(A(mask));
% 
%                 %%%Calculate centersX of max receptive field
% 
%                 %%%Calcutate centersY
% 
%                 %%%Calculate locality index           
% 
%                 %index(i,ui) = maxMean/max(meanA,[],'all') - outsideMean/max(meanA,[],'all');
% 
%                 index(3,ui) = maxMean;
%                 index(4,ui) = outsideMean;
% 
%                 if maxMean>maxAllEM
%                     maxAllEM = maxMean;
%                 end
% 
% 
%                 indexFinal(3,ui) = d;
%                 indexFinal(4,ui) = respU(u);
%                 indexFinal(5,ui) = in;
% 
%                 ui = ui+1;
%            
%             end
% 
%             %NEM
%             indexFinal(1,ui-4:ui-1) = index(1,ui-4:ui-1)./maxAllNEM - index(2,ui-4:ui-1)./maxAllNEM;
% 
%             %EM
%             indexFinal(2,ui-4:ui-1) = index(3,ui-4:ui-1)./maxAllEM - index(4,ui-4:ui-1)./maxAllEM;
% 
% 
%     end
% 
%     in = in+1;
% 
%     indexSummary = [indexSummary indexFinal]; 


if perDirection ==0
    %%%% Load real responses:
    %%noEyeMoves
    NMRF = squeeze(load(sprintf('NEM-RFuST-Q1-Div-X-%s',NP.recordingName)).RFuST);
    RFrectGrid =  load(sprintf('RFuStatic-%s',NP.recordingName)).RFu;
    %%EyeMoves
    RF = squeeze(load(sprintf('RFuST-Q1-Div-X-%s',NP.recordingName)).RFuST);
    RFReal = {NMRF,RF};

    %%% Load suffled responses:
    %%NoEyeMoves
    NMRFs = squeeze(load(sprintf('NEM-RFuShuffST-Q1-Div-X-%s',NP.recordingName)).RFuShuffST);
    %RFrectGridShuff = load(sprintf('RFuStaticShuff-%s',NP.recordingName)).RFu;
    %%EyeMoves
    RFs = squeeze(load(sprintf('RFuShuffST-Q1-Div-X-%s',NP.recordingName)).RFuShuffST);
    RFShuffle = {NMRFs,RFs};


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate tuning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% indexes

    %1. Bin the receptive field into offsetxoffset grid and
    %take the mean within each grid.

    perc = 90; %Tope percentile selected to calculate tuning

    rowT = cell(1,numel(respU));
    colT = cell(1,numel(respU));

    SpatialTuningIndex = zeros(4,numel(respU));
    %
    % figure;tiledlayout(9,10,"TileSpacing","tight")
    %
    % MrC = ConvBurstMatrix(Mr,fspecial('gaussian',[1 10],3),'same');
    %shuffC = ConvBurstMatrix(squeeze(mean(shuffledData,4)),fspecial('gaussian',[1 10],3),'same');

    for em =1:2

        RFuST = RFReal{em};

        RFuShuffST = RFShuffle{em};

        meanY = zeros(1,length(respU));
        meanX = zeros(1,length(respU));

        for u =1:length(respU)

            % Example 54x54 matrix
            rfi = squeeze(RFuST(:,:,u));

            % Define block size
            blockSize = redCoorY/offsetN;

            % Reshape and compute the mean for each 9x9 block
            meanGrid = blockproc(rfi, [blockSize blockSize], @(x) mean(x.data(:)));

            % Find mean of top blocks and rest, and the location of
            % the top blocks
            TOPper = prctile(meanGrid(:), perc);

            Rtop = mean(meanGrid(meanGrid >= TOPper));

            [rowT, colT] = find(meanGrid >= TOPper);

            meanY(u) = mean(rowT)/size(meanGrid,1);
            meanX(u) = mean(colT)/size(meanGrid,1);

            Rrest = mean(meanGrid(meanGrid < TOPper));

            if isnan(Rrest) %% In case top percentile is values bigger than 0 (sparse unit)
                Rrest = 0;
            end

            % Calculate first term

            realTerm = (Rtop - Rrest)/mean(meanGrid,'all');

            %Do the same for the shuffled

            shuffTerm =zeros(1,numel(respU));

            for s = 1:nShuffle

                rfiShuff = squeeze(RFuShuffST(:,:,u,s));

                meanGridShuff = blockproc(rfiShuff, [blockSize blockSize], @(x) mean(x.data(:)));
                %
                %                         % Use the same positions of the top blocks in the
                %                         % shuffled data
                %
                %                         RtopS = mean(meanGridShuff(rowT{u},colT{u}),'all');
                %
                %                         Mask = true(size(meanGridShuff));
                %
                %                         Mask(rowT{u},colT{u}) = false;
                %
                %                         RrestS = mean(meanGridShuff(Mask));

                % Find mean of top blocks and rest, and the location of
                % the top blocks
                TOPper = prctile(meanGridShuff(:), perc);

                RtopS = mean(meanGridShuff(meanGridShuff >= TOPper));

                RrestS = mean(meanGridShuff(meanGridShuff < TOPper));

                if isnan(Rrest)
                    Rrest = 0;
                end

                % Calculate first term

                shuffTerm(s) = (RtopS - RrestS)/mean(meanGridShuff,'all');

            end

            SpatialTuningIndex(em,u) = realTerm-median(shuffTerm);

            SpatialTuningIndex(3,u) = respU(u);

            SpatialTuningIndex(4,u) = in;

            %     nexttile
            %     imagesc(rfi);title(string(SpatialTuningIndex(u)))
            %     axis off;

            %imagesc(squeeze(shuffC(:,u,:)));title(string(shuffTerm(u)))

            %                     imagesc(squeeze(MrC(:,respU(u),:)));title(string(SpatialTuningIndex(u)))
            %                     yline([10:10:size(Mr,1)],'w','LineWidth',1);
            %                     yline([90:90:size(Mr,1)],'w','LineWidth',3);

        end %%% end units

    end %%% end eye moves or not

else %%%Do per direction

     %%%% Load real responses:
    %%noEyeMoves
    NMRF = squeeze(load(sprintf('NEM-RFuSTDirSize-Q1-Div-X-%s',NP.recordingName)).RFuSTDirSize);
    %%EyeMoves
    RF = squeeze(load(sprintf('RFuSTDirSize-Q1-Div-X-%s',NP.recordingName)).RFuSTDirSize);

    RFReal = {NMRF,RF};

    %%% Load suffled responses:
    %%NoEyeMoves
    NMRFs = squeeze(load(sprintf('NEM-RFuShuffST-Q1-Div-X-%s',NP.recordingName)).RFuShuffST);
    %%EyeMoves
    RFs = squeeze(load(sprintf('RFuShuffST-Q1-Div-X-%s',NP.recordingName)).RFuShuffST);

    RFShuffle = {NMRFs,RFs};

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Calculate tuning
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% indexes

    %1. Bin the receptive field into offsetxoffset grid and
    %take the mean within each grid.

    perc = 90; %Tope percentile selected to calculate tuning

    rowT = cell(1,numel(respU));
    colT = cell(1,numel(respU));

    SpatialTuningIndex = zeros(5,numel(respU)*4);
    %
    % figure;tiledlayout(9,10,"TileSpacing","tight")
    %
    % MrC = ConvBurstMatrix(Mr,fspecial('gaussian',[1 10],3),'same');
    %shuffC = ConvBurstMatrix(squeeze(mean(shuffledData,4)),fspecial('gaussian',[1 10],3),'same');

    for em =1:2

        ui =1;

        RFuST = RFReal{em};

        RFuShuffST = RFShuffle{em};

        for u =1:length(respU)

            for d =1:4

                % Example 54x54 matrix
                rfi = squeeze(RFuST(d,:,:,u));

                % Define block size
                blockSize = redCoorY/offsetN;

                % Reshape and compute the mean for each 9x9 block
                meanGrid = blockproc(rfi, [blockSize blockSize], @(x) mean(x.data(:)));

                % Find mean of top blocks and rest, and the location of
                % the top blocks
                TOPper = prctile(meanGrid(:), perc);

                Rtop = mean(meanGrid(meanGrid >= TOPper));

                [rowT{u}, colT{u}] = find(meanGrid >= TOPper);

                Rrest = mean(meanGrid(meanGrid < TOPper));

                if isnan(Rrest) %% In case top percentile is values bigger than 0 (sparse unit)
                    Rrest = 0;
                end

                % Calculate first term

                realTerm = (Rtop - Rrest)/mean(meanGrid,'all');

                if d == 1 %%Calculate shuffled indez once per direction

                    %Do the same for the shuffled

                    shuffTerm =zeros(1,numel(respU));

                    for s = 1:nShuffle

                        rfiShuff = squeeze(RFuShuffST(:,:,u,s));

                        meanGridShuff = blockproc(rfiShuff, [blockSize blockSize], @(x) mean(x.data(:)));

                        % Find mean of top blocks and rest, and the location of
                        % the top blocks
                        TOPper = prctile(meanGridShuff(:), perc);

                        RtopS = mean(meanGridShuff(meanGridShuff >= TOPper));

                        RrestS = mean(meanGridShuff(meanGridShuff < TOPper));

                        if isnan(Rrest)
                            Rrest = 0;
                        end

                        % Calculate first term

                        shuffTerm(s) = (RtopS - RrestS)/mean(meanGridShuff,'all');

                    end

                end

                SpatialTuningIndex(em,ui) = realTerm-median(shuffTerm);

                SpatialTuningIndex(3,ui) = respU(u);

                SpatialTuningIndex(4,ui) = in;

                SpatialTuningIndex(5,ui) = d;

                ui = ui+1;

            end %%% end directions

        end %%% end units

    end %%% end eye moves or not

end

indexSummary = [indexSummary SpatialTuningIndex];

locationRF = [locationRF [meanX;meanY]];

in = in+1;

end %%% end insertions



%% Colors per insertion

% Number of colors
nColors = in-1;

% Generate hues evenly spaced around the color wheel
H = linspace(0, 1, nColors+1); 
H = H(1:end-1); % Remove last to avoid duplication at 1 and 0

% Keep saturation low for pale colors
S = 0.5 * ones(1, nColors); 

% Keep value high for brightness
V = 0.7 * ones(1, nColors); 

% Convert HSV to RGB
colorMatrix = hsv2rgb([H' S' V']);

% Display the colors
figure;
colormap(colorMatrix);
c = colorbar;
caxis([0 1])
c.Ticks = linspace(1/(in-1),1,in-1)-(1/(in-1))/2; 
c.TickLabels = 1:in-1;
c.Limits
title('Pale Differentiated Colors');
legend()


colors = colorMatrix(indexSummary(4,:),:);
%% Scatter

figure;scatter(locationRF(1,:),locationRF(2,:),7,colors,"filled",'MarkerFaceAlpha',0.7)
set(gca, 'YDir', 'reverse');
cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
title('Locality index')
axis equal
print(gcf, 'RFsLocations.pdf', '-dpdf', '-r300', '-vector');


%% Histogram


figure;histogram(indexSummary(2,:),'FaceColor','k','FaceAlpha',0.5)
hold on;
histogram(indexSummary(1,:),'FaceColor','b','FaceAlpha',0.5)
ylabel('Neuron count')
xlabel('L')
cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
print(gcf, 'LocalityIndexDistribution.pdf', '-dpdf', '-r300', '-vector');
cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
title('Locality index')
print(gcf, 'LocalityIndex.pdf', '-dpdf', '-r300', '-vector');

%% Plot NEM vs EM %Color



fig = figure;
scatter(indexSummary(1,:),indexSummary(2,:),7, colors,"filled",'MarkerFaceAlpha',0.4);
xlabel('C.Eye Mov')
ylabel('NC. Eye Mov')
hold on
plot([0,max(indexSummary(1:2,:),[],'all')],[0,max(indexSummary(1:2,:),[],'all')],'LineWidth',1,'Color','k')
axis equal
fig.Position = [1141         261         214         140];
xlim([0,max(indexSummary(1:2,:),[],'all')])
ylim([0,max(indexSummary(1:2,:),[],'all')])
cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
title('Locality index')
print(fig, 'LocalityIndex.pdf', '-dpdf', '-r300', '-vector');


%% Plot swarm per directions


%%Create a round histogram divided into 4 direction. Within each direction,
%%DSI values are sorted

angles = indexSummary(5,:);
L = indexSummary(1,:);

fig=figure;swarmchart(categorical(angles),(L),10, colors,'filled','MarkerFaceAlpha',0.8);

set(gcf,'Color','w');%
ax = gca; % Get current axis
ax.XAxis.FontSize = 7; % Set font size of x-axis tick labels
ylabel('log(L)')
yline(0,'LineWidth',1,'Color','k')
set(gca, 'YScale', 'log')

%fig.Position = [1269         521         189         146];

print(gcf, 'swarmPlotLocalIndexPerAngle.pdf', '-dpdf', '-r300', '-vector');

%% Plot histogram of neuron counts per direction  that have L>0.5
fig = figure;
Dirs= {};
for i =1:in-1
Dirs{i} = categorical(indexSummary(3,indexSummary(1,:)>0.5 & indexSummary(5,:)==i));

end

figure;histogram(Dirs);


bar([1:4], tuningCurve(u,:).*1000, 'FaceColor', 'k','FaceAlpha',0.5); % Bar plot for spike rates

hold on;

fig.Position = [  997   266   149   135];

figure;histogram


%% Plot several examples PV97 & PV35
cd('D:\Mark_S13\Desktop\receptiveFieldExamples');
file = dir ('D:\Mark_S13\Desktop\receptiveFieldExamples');
filenames = {file.name};

append_pdfs(sprintf('RF-%s',NP.recordingName),filenames{:})

figure;plot(indexSummary(1,:),indexSummary(2,:))



% crossCor 
% 
% cc = xcorr2(offsetTemplate,template);
% [max_cc, imax] = max(abs(cc(:)));
% [ypeak, xpeak] = ind2sub(size(cc),imax(1));
% corr_offset = [(ypeak-size(template,1)) (xpeak-size(template,2))];