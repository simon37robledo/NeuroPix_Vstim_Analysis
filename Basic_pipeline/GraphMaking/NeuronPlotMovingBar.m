%%%%% Neuron plot control

%%% Optional inputs:
%
% 'SelectRand' = integer. Should be less than total trials. Select random trials to signal in the raster, plot the raw data, and eye movements.
% 'SelectRand' = integer. Should be less than total trials. Select random trials to signal in the raster, plot the raw data, and eye movements.
% 'SelectRand' = integer. Should be less than total trials. Select random trials to signal in the raster, plot the raw data, and eye movements.
% 'SelectRand' = integer. Should be less than total trials. Select random trials to signal in the raster, plot the raw data, and eye movements.
%%
function NeuronPlotMovingBar(data,ex,eNeuron,varargin)

% Create parser object
p = inputParser;

%% Add required inputs
addRequired(p, 'data');
addRequired(p, 'ex');
addRequired(p, 'eNeuron');

%% Optional Name-Value pairs with defaults

%Raster
addParameter(p, 'Raster', 0);
addParameter(p, 'TrialNumber', 1);
addParameter(p, 'SelectRand', 0);
addParameter(p, 'window', 300); %window (ms)
addParameter(p, 'start', -50); %Start respective to stimulus onset (ms)

%Tuning curve
addParameter(p, 'tuningPlot',0)

%Receptive field
addParameter(p, 'ReField', 0);
addParameter(p, 'oneDir', 0);
addParameter(p, 'EyeMovements', 0);
addParameter(p, 'noEyeMoves', 0);
addParameter(p, 'DivisionType', 'Mode');

%Save plot
addParameter(p, 'savePlot', 0);
addParameter(p, 'saveDir', '\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure');

%% Parse inputs
parse(p, data,ex,eNeuron, varargin{:});

% Extract values

%Raster
Raster = p.Results.Raster;
TrialNumber = p.Results.TrialNumber;
SelectRand = p.Results.SelectRand;
window = p.Results.window;
start = p.Results.start;

%Tuning curve
tuningPlot = p.Results.tuningPlot;

%Receptive field
ReField = p.Results.ReField;
oneDir = p.Results.oneDir;
EyeMovements = p.Results.EyeMovements;
noEyeMoves = p.Results.noEyeMoves;
DivisionType  =p.Results.DivisionType;

%save plot
savePlot = p.Results.savePlot;
saveDir = p.Results.saveDir;


% Loop through varargin as Name-Value pairs
for i = 1:2:length(varargin)
    name = varargin{i};   % Extract name
    value = varargin{i+1}; % Extract value
    fprintf("x%s = %s\n", name, mat2str(value)); % Print output
end

%% Load NP class %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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


%% %%%%%%%%%%%%%%% load stimulus properties %%%%%%%%%%%%%%%%%%%%%%%%%%%



cd(NP.recordingDir)

%2. Extract moving ball statistics
patternIndex = strfind(string(NP.recordingDir), "\catgt");

endIndex = patternIndex(1)-1;
stimDir = string(NP.recordingDir);
stimDir = extractBetween(stimDir,1,endIndex);

file = dir (stimDir);
filenames = {file.name};

file = dir (stimDir);
filenames = {file.name};
ballFiles = filenames(contains(filenames,"linearlyMovingBar"));


if isempty(ballFiles)
    %disp()
    w= sprintf('No moving bar files where found in %s. Skipping into next experiment.',NP.recordingName);
    warning(w)

end

    directions = [];
    lengthBar = [];
    widthBar = [];
    speeds = [];

    j =1;

    if size(ballFiles) ~= [0 0]

        for i = ballFiles %Extract stim properties
            ball= load(stimDir+"\"+string(i));

            directions = [directions cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'directions'))))];

            widthBar = [widthBar cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'barWidth'))))];

            lengthBar = [lengthBar cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'barLength'))))];

            speeds = [speeds cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'speeds'))))];

            interStimStats = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

%% %%%%%%%%%%%%%%%%%%%% load neuron properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = NP.convertPhySorting2tIc(NP.recordingDir);
label = string(p.label');
goodU = p.ic(:,label == 'good');
phy_IDg = p.phy_ID(label == 'good');

% try
%     NeuronVals = load(sprintf('OneSpeed-%d-NeuronRespCat-%s',max(speeds),NP.recordingName)).NeuronVals;
% catch
%     NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;
% end

%% %%%%%%%%%%%%%% Stim specific parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load Triggers (diode)
Ordered_stims= strsplit(data.VS_ordered{ex},',');
%containsMB = cellfun(@(x) contains(x,'MB'),Ordered_stims);
containsMB = strcmp(Ordered_stims, 'MBR');
ttlInd = find(containsMB);
newDiode =0;
[stimOn stimOff onSync offSync] = NPdiodeExtract(NP,newDiode,1,"MBR",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
[stimOn stimOff onSync offSync] = NPdiodeExtract(NP,0,1,"MBR",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A


%When dealing with different speeds, save different stim durs, and create M for each speed
A = [stimOn directions' speeds'];
[C indexS] = sortrows(A,[2 3]);

B = [stimOff directions' speeds'];
[Coff indexSo] = sortrows(B,[2 3]);

stimInter= mean(stimOn(2:end)-stimOff(1:end-1));

if includeOnespeed

    A =  A(A(:,3) == max(A(:,3)),:);

    B =  B(B(:,3) == max(B(:,3)),:);

    [C indexS] = sortrows(A,[2 3]);

    [Coff indexSo] = sortrows(B,[2 3]);

    stimOn = A(:,1);
    stimOff = B(:,1);

    %digon  =onDigital(A(:,5) == max(A(:,5)));

    stimInter = interStimStats;

end


uDir = unique(directions);
direcN = length(uDir);
uSpeed = unique(speeds);
speedN = length(uSpeed);
%offsetN = length(unique(offsets));
nT = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials'))));
trialDivision = nT/(length(uDir)*length(uSpeed)); %Number of trials per unique conditions


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%% Build matrix %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
directimesSorted = C(:,1)';

stimDur = mean(-stimOn+stimOff);

preBase = round(stimInter-stimInter/4);

bin = 1;
bin2 =15;

[Mr2] = BuildBurstMatrix(goodU(:,eNeuron),round(p.t/bin2),round((directimesSorted-preBase)/bin2),round((stimDur+preBase*2)/bin2));

[nT,nN,nB] = size(Mr2);


N_bootstrap = 1000;

pvals = load(sprintf('pvalsBaselineBootMBR-1000-%s',N_bootstrap,NP.recordingName),'pvalsResponse').pvalsResponse;

pvals = [eNeuron;pvals(eNeuron)]; %%%Why unit 59 has a pvalue of 0. 


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot raster MB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ur = 1;

for u = eNeuron

    


    if Raster

        fig =  figure;

        title(sprintf('U.%d-Unit-phy-%d-p-%d',u,phy_IDg(u),pvals(2)));


        subplot(20,1,[7 18]);
        imagesc(squeeze(Mr2).*(1000/bin2));colormap(flipud(gray(64)));
        %Plot stim start:
        xline(preBase/bin2,'k', LineWidth=1.5)
        %Plot stim end:
        xline(stimDur/bin2+preBase/bin2,'k',LineWidth=1.5)


        caxis([0 1])
%         dirStart = C(1,2);
%         offStart = C(1,3);
%         for t = 1:nT*mergeTrials
%             if dirStart ~= C(t,2)
%                 yline(t/mergeTrials-0.5,'k',LineWidth=1.5);
%                 dirStart = C(t,2);
%             end
%             if offStart ~= C(t,3)
%                 yline(t/mergeTrials-0.5,'k',LineWidth=0.5);
%                 offStart = C(t,3);
%             end
% 
%         end

        hold on


        xticklabels([])
        xlim([0 round(stimDur+preBase*2)/bin2])
        xticks([0 preBase/bin2:300/bin2:(stimDur+preBase*2)/bin2 (round((stimDur+preBase*2)/100)*100)/bin2])
        xticklabels([]);

        yt =[0];
        for d = 1:direcN
            yt = [yt [1:trialDivision*2*sizeN:(nT/direcN)-1+trialDivision*sizeN]+max(yt)+trialDivision-1];

        end

        yt = yt(2:end);
        yticks(yt)

        yticklabels(repmat([trialDivision:trialDivision*2*sizeN:(nT/direcN)-1+trialDivision*sizeN],1,direcN))

        ax = gca; % Get current axes
        ax.YAxis.FontSize = 7; % Change font size of y-axis tick labels


        %%%%Select 10 random trials where Z-score is above the average

        SelectRand = 0;

        if SelectRand > 0

            ZscoreRaster = load(sprintf('ZscoreRaster-%d-%s',N_bootstrap,NP.recordingName)).ZscoreRaster;

            mZscoreRaster = mean(ZscoreRaster(:,u,:),3);

            trialsHigherSpiking = find(mZscoreRaster>mean(mZscoreRaster)+2*std(mZscoreRaster));

            trials = trialsHigherSpiking(sort(randperm(length(trialsHigherSpiking),10)));

            %%Hightlight trials in raster:

            colors = turbo(10);
            xlimits = xlim;

            yl = yline(trials,'LineWidth',2.5);

            for i = 1:numel(trials)

                yl(i).Color = [colors(i,:) 0.2];

            end

        else %%Select highest window stim type

            j =1;
            meanMr = zeros(1,nT/trialDivision);
            for i = 1:trialDivision:nT
               
                meanMr(j) = mean(Mr2(i:i+trialDivision-1,:),'all');

                j = j+1;
            end


            %[maxResp,maxRespIn]= max(NeuronVals(u,:,1));

            [maxResp,maxRespIn]= max(meanMr);

            maxRespIn = maxRespIn-1;

            trials = maxRespIn*trialDivision+1:maxRespIn*trialDivision + trialDivision;

            colors = turbo(trialDivision);

            y1 = maxRespIn*trialDivision + trialDivision;
            y2 = maxRespIn*trialDivision;

            %Figure paper
%             start = -50;
%             window = 300;

            patch([(preBase+start)/bin2 (preBase+start+window)/bin2 (preBase+start+window)/bin2 (preBase+start)/bin2],...
                [y2 y2 y1 y1],...
                'b','FaceAlpha',0.3,'EdgeColor','none')

             patch([1 (preBase*2+stimDur)/bin2 (preBase*2+stimDur)/bin2 1],...
                [y2 y2 y1 y1],...
                'b','FaceAlpha',0.2,'EdgeColor','b')


        end

        %TrialNumber = 8;
        %Mark selected trial
        RasterTrials = trials(TrialNumber);
% 
%         spikes1 = Mr2(RasterTrials,:);
%         
%         hold on
%         %plot(zeros(1,numel(RasterTrials)),RasterTrials,'*')
%        
%         spikeLoc = find(spikes1 >0);
%         plot(spikeLoc, repmat(RasterTrials,numel(spikeLoc),1),'.','MarkerSize',10,'MarkerFaceColor','r','MarkerEdgeColor','r')

        %%%%%% Plot PSTH
        %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        subplot(20,1,[19 20])

        MRhist = BuildBurstMatrix(goodU(:,u),round(p.t),round((directimesSorted-preBase)),round((stimDur+preBase*2)));

        MRhist = squeeze(MRhist(trials,:,:));

        [nT2,nB2]= size(MRhist);

        spikeTimes = repmat([1:nB2],nT2,1);

        spikeTimes = spikeTimes(logical(MRhist));
        % Define bin edges (adjust for resolution)
        binWidth = 125; % 10 ms bins
        edges = [1:binWidth:round((stimDur+preBase*2))]; % Adjust time window as needed

        % Compute histogram
        psthCounts = histcounts(spikeTimes, edges);

        % Convert to firing rate (normalize by bin width)
        psthRate = (psthCounts / (binWidth * nT2))*1000;

        b=bar(edges(1:end-1), psthRate,'histc');
        b.FaceColor = 'b';
        b.FaceAlpha = 0.3;
        b.MarkerEdgeColor = "none";
        ylabel('Spikes/sec','FontSize',10);
        xlabel('Seconds','FontSize',10);
        xlim([0 round((stimDur+preBase*2)/100)*100])
        ylim([0 max(psthRate)+std(psthRate)])

        xticks([0 preBase:300:(stimDur+preBase*2) round((stimDur+preBase*2)/100)*100])

        xline([preBase stimDur+preBase],'LineWidth',1.5)

        xticklabels([-(preBase) 0:300:round((stimDur/100))*100 round((stimDur/100))*100 + preBase]./1000)
      
        ax = gca; % Get current axes
        ax.YAxis.FontSize = 7; % Change font size of y-axis tick labels
        ax.XAxis.FontSize = 7;


        %set(gca, 'XLim', [-0.5, 1]); % Match raster plot x-axis


        %%%%PLot raw data several trials one
        %%%%channel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            chan = goodU(1,u);
            
            subplot(20,1,[1 3])
           
            startTimes = directimesSorted(RasterTrials)+start;

            freq = "AP"; %or "LFP"

            type = "line"; %or heatmap

            spikes = squeeze(BuildBurstMatrix(goodU(:,u),round(p.t),round(startTimes),round((window))));

            [fig, mx, mn] = PlotRawData(fig,NP,chan,startTimes,window,freq,type,preBase,spikes');

            xlabel(string(chan))
            xline(-start/1000,'LineWidth',1.5)
            xticklabels([])
            xlabel([]);xticks([])
            ylabel('uV')
                
        
        %%%%%%%%%%% Plot raster of selected trials
        %%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%
            subplot(20,1,[4 5])

            %RasterTrials =  trials(2); %Trial selected above

            RasterTrials =  trials;

            bin3 = 2;
            trialM = BuildBurstMatrix(goodU(:,eNeuron),round(p.t/bin3),round((directimesSorted+start)/bin3),round((window)/bin3));
            
            if numel(RasterTrials) == 1         
                TrialM = squeeze(trialM(RasterTrials,:,:))';
            else               
                TrialM = squeeze(trialM(trials,:,:));   
            end
            

            TrialM(TrialM~=0) = 0.3;
            spikes1 = TrialM(TrialNumber,:);
            spikeLoc = find(spikes1 >0);
            TrialM(TrialNumber,spikeLoc) = 1;

            %Select offset in which selected trial belongs (10 trials)
            imagesc(TrialM);colormap(flipud(gray(64)));
            %caxis([0 1])
            xline([preBase stimDur+preBase],'LineWidth',1.5)
            ylabel([sprintf('%d trials',numel(trials))])
            xticks([1 abs(start)/bin3:20/bin3:window/bin3])
            xticklabels([start 0:20:window-abs(start)])
            xlabel('Milliseconds')
            set(gca,'FontSize',7)

            xline(-start/bin3,'LineWidth',1.5)
            
            xline(spikeLoc,'LineWidth',1,'Color','b','Alpha',0.3)

            yline(TrialNumber,'LineWidth',3,'Color','b','Alpha',0.3) %Mark trial
%   
            fig.Position = [680     5   296   9734];
% 
%         limsY = ylim;
%         ylim([limsY(1)-mx limsY(2)+mx]); %%%Set limits using max and min spike values

        if savePlot
            cd(saveDir)
            print(gcf, sprintf('%s-MovBall-SelectedTrials-eNeuron-%d.pdf',NP.recordingName,u), '-dpdf', '-r300', '-vector');
        end

        %if max(RasterTrials) > 0 %%trials organized by direction


        %end

    end
        
    



    %% %%%%%%%%%Plot receptive field per direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if ReField

        close all

        cd(NP.recordingDir)

        %Parameters
        %check receptive field neurons first
        respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;

        respU = find(respNeuronsMB<0.05);

        ru = find(respU == u);
        %
        reduceFactor =20;
        coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));
        reduceFactor = min([20 min(sizes)]); %has to be bigger than the smallest ball size
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

        if oneDir

            if noEyeMoves %%%mode quadrant
                RFu = squeeze(load(sprintf('NEM-RFuST-Q1-Div-X-%s',NP.recordingName)).RFuST);
            else
                RFu = squeeze(load(sprintf('RFuST-Q1-Div-X-%s',NP.recordingName)).RFuST);
            end

            %%%Filter with gaussian:

            TwoDGaussian = fspecial('gaussian',floor(size(RFu,2)/(9/2)),redCoorY/offsetN); %increase size of gaussian by 100%.
            %RFuSTDirFilt = zeros(size(RFuSTDir));

            RFuFilt= zeros(size(RFu));

            for ui =1:size(RFu,3) %units

                slice = squeeze(RFu(:,:,ui));

                slicek = conv2(slice,TwoDGaussian,'same');

                RFuFilt(:,:,ui) =slicek;
            end

            figRF=figure;
            imagesc((squeeze(RFuFilt(:,:,ru))));

          caxis

            c = colorbar;
            title(c,'spk/s')

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

            if savePlot
                cd(saveDir)
                figRF.Position = [ 680   577   156   139];
                if noEyeMoves
                    print(figRF, sprintf('%s-NEM-MovBall-ReceptiveField-eNeuron-%d.pdf',NP.recordingName,u), '-dpdf', '-r300', '-vector');
                else
                    print(figRF, sprintf('%s-MovBall-ReceptiveField-eNeuron-%d.pdf',NP.recordingName,u), '-dpdf', '-r300', '-vector');
                end
            end

        else %%%% Plot receptive field per direction

        if noEyeMoves

            if string(DivisionType) == "XY"

                file = dir (NP.recordingDir);
                filenames = {file.name};
                cd(NP.recordingDir)

                names = {'X','Y'};

                files = filenames(contains(filenames,"NEM-RFuSTDirSizeFilt")& contains(filenames,"Div"));


                for n = 1:2

                     filesD = filenames(contains(files,names{n}));

                    for f = 1:numel(filesD)
                        rf{1,f} = load(sprintf('NEM-RFuSTDirSizeFilt-Q%d-Div-%s-%s',f,names{n},NP.recordingName)).RFuSTDirSizeFilt;

                        rf{2,f} = sprintf('NEM-RFuSTDirSizeFilt-Q%d-Div-%s-%s',f,names{n},NP.recordingName);
                    end
  
                    RF{1,n} = rf;

                    RF{2,n} = sprintf('Division along %s',names{n});
                end

                %RF contains two cells, one for each type of eye movement
                %division (along X or along Y)

%                 q=1;
%                 RF{1} = load(sprintf('NEM-RFuSTDirSizeFilt-Q%d-%s',q,NP.recordingName)).RFuSTDirSizeFilt; %Size and dir
%                 q=2;
%                 RF{2} = load(sprintf('NEM-RFuSTDirSizeFilt-Q%d-%s',q,NP.recordingName)).RFuSTDirSizeFilt; %Size and dir

            else

                rf{1} = load(sprintf('NEM-RFuSTDirSizeFilt-Q1-Div-X-%s',NP.recordingName)).RFuSTDirSizeFilt; %Size and dir

                RF{1} = rf;

                names{1} ="";
            end
        else
            q=1;
            rf{1} = load(sprintf('RFuSTDirSizeFilt-%s',NP.recordingName)).RFuSTDirSizeFilt; %Size and dir

            RF{1} = rf;

            names{1} ="";
        end   

        %%%% find max and min of colorbar limits

        cMax = -inf;
        cMin = inf;
        for n = 1:size(RF,2)

            rf = RF{1,n};
            for q = 1:size(rf,2)
                RFuSTDirSizeFilt = rf{1,q};

                if size(RFuSTDirSizeFilt,2) >1 %If there is more than on size %%Solve for 2 quadrants
                    %%Select size closest to 120
                    [minS indx] = min(abs(unique(sizes)-120));
                    RFuSTDirFilt = squeeze(RFuSTDirSizeFilt(:,indx,:,:,:));
                else
                    RFuSTDirFilt = squeeze(RFuSTDirSizeFilt);
                end

                for d = 1:direcN

                    if cMax < max(RFuSTDirFilt(:,:,:,ru),[],'all')

                        cMax = max(RFuSTDirFilt(:,:,:,ru),[],'all');
                    
                    end

                      if cMin > min(RFuSTDirFilt(:,:,:,ru),[],'all')

                        cMin = min(RFuSTDirFilt(:,:,:,ru),[],'all');                   
                      
                      end

                end
            end

        end

         % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figRF = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full screen figure;
        NeuronLayout = tiledlayout(size(RF,2),size(RF,2),"TileSpacing","tight","Padding","tight");

        tl = 1;
        for n = 1:size(RF,2)

            rf = RF{1,n};

            for q = 1:size(rf,2)

                if names{n} == 'X'
                    qNames = {'Right','Left'};

                elseif names{n} == 'Y'
                    qNames = {'Down','Up'};
                else
                    qNames = {"-"};
                end

                RFuSTDirSizeFilt = rf{1,q};

                if size(RFuSTDirSizeFilt,2) >1 %If there is more than on size %%Solve for 2 quadrants
                    %%Select size closest to 120
                    [minS indx] = min(abs(unique(sizes)-120));
                    RFuSTDirFilt = squeeze(RFuSTDirSizeFilt(:,indx,:,:,:));
                else
                    RFuSTDirFilt = squeeze(RFuSTDirSizeFilt);
                end


                xqLayout=tiledlayout(NeuronLayout,direcN/2,direcN/2,"TileSpacing","compact");
                xqLayout.Layout.Tile = tl;

                %DirLayout=tiledlayout(direcN/2,direcN/2,"TileSpacing","tight","Padding","none");


                for d = 1:direcN

                    nexttile(xqLayout);

                    if d ==1 || d==3
                        %imagesc(rot90(squeeze(RFuSTDirFilt(d,:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY,ru)),2));

                        imagesc((squeeze(RFuSTDirFilt(d,:,:,ru))));
                        
%                         c = colorbar;
                        2+2
                    else
                        imagesc((squeeze(RFuSTDirFilt(d,:,:,ru))));
                        
                    end

                    if d==2
                    c = colorbar;
                    title(c,'Spike Rate (Hz)')
                    end

                    

                    caxis([cMin cMax]);
                    
                    colormap('turbo')
                    title(string(uDir(d)))
                    
                    %xlim([(redCoorX-redCoorY)/2 (redCoorX-redCoorY)/2+redCoorY])
                    xt = xticks;%(linspace((redCoorX-redCoorY)/2,(redCoorX-redCoorY)/2+redCoorY,offsetN*2));
                    xt = xt((1:2:numel(xt)));
                    xticks(xt);
                    xticklabels(round(theta_x(1,xt)))
                    
                    yt = yticks;
                    yt = yt(1:2:numel(yt));
                    yticks(yt);
                    yticklabels(round(theta_y(yt,1)))
                    %                     xlabel('X degrees')
                    %                     ylabel('Y degrees')
                    if d == 1 || d==2

                        xticks([])

                    end

                    if d==2 || d==4

                        yticks([])

                    end

%%%%Use same colorbar scale

                    axis equal tight


                end

                ti = title(xqLayout,sprintf('nUnit-%d-Div-%s-Q-%s\n\',respU(ru),names{n},qNames{q}),'FontSize',10);

                %ti.Position(2) = ti.Position(2) + 0.05;

                tl = tl+1;

                %         %%%%%%%%% plot rasters used in noEyeMoves RFs
                %
                %         spikeSum = spikeSums{q};
                %
                %         eNspkS = squeeze(spikeSum(:,ru,:));
                %         figure;imagesc(eNspkS);colormap(flipud(gray(64)));
                %         %     rowsWithNaN = find(any(isnan(eNspkS), 2));
                %         %     yline(rowsWithNaN,'g')
                %         yline(trialDivision*sizeN:trialDivision*sizeN:size(spikeSums{q},1));
                %         yline(trialDivision*offsetN*sizeN:trialDivision*offsetN*sizeN:size(spikeSums{q},1)-1,'r','LineWidth',5)
                %         title(string(q))

            end %%%%%%%%%%%%%%%%

           

        end

        if savePlot
            cd(saveDir)
            figRF.Position = [  0.2328125                     0.315                0.23515625                   0.38125];
            if noEyeMoves
                print(gcf, sprintf('%s-NEM-MovBall-ReceptiveField-Div-%s-Q-%d-eNeuron-%d.pdf',NP.recordingName,names{n},q,u), '-dpdf', '-r300', '-vector');
            else
                print(gcf, sprintf('%s-MovBall-ReceptiveField-Div-%s-Q-%d-eNeuron-%d.pdf',NP.recordingName,names{n},q,u), '-dpdf', '-r300', '-vector');
            end
        end
        
        end %%%End onDir


    end

    if EyeMovements


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%Plot eye movements colorcoded for these selected trials:

        %         %%%%%% Construct stimType matrix for eye movement plot.
        %         stimType = zeros(length(A),5); %3 plus number of example neuron
        %         stimType(:,1) = A(:,1);
        %         stimType(:,2) = A(:,1)+stimDur;
        %         stimType(:,3) = A(:,2);
        %         stimType(:,4) = A(:,3);
        %
        %         %Get response strenght of specific neurons and save it in stimType
        %         [MrNoSort] = BuildBurstMatrix(goodU(:,u),round(p.t/bin),round((stimOn'-preBase)/bin),round((stimDur+preBase*2)/bin)); %response matrix
        %         ResponseStrength= mean(MrNoSort(:,u,round(preBase/bin):round((preBase+stimDur)/bin)),3); %For PV35_3
        %         stimType(:,end) = ResponseStrength;

        if SelectRand
            timesnips = [directimesSorted(trials);directimesSorted(trials)+stimDur];

            EyePositionAnalysis(NP,data.Eye_video_dir{ex},21,1,0,0,timesnips,1,0)
        elseif DivisionType == "XY" %%%divided

            EyePositionAnalysis(NP,data.Eye_video_dir{ex},21,1,0,0,0,1,1)

        else %% mode type of division

            EyePositionAnalysis(NP,data.Eye_video_dir{ex},31,1,0,0,0,1,0)

        end

        figure

        if savePlot
            cd(saveDir)
            print(gcf, sprintf('%s-MovBall-EyeMovs-SelectedTrials-eNeuron-%d.pdf',NP.recordingName,u), '-dpdf', '-r300', '-vector');
        end

    end

    if tuningPlot


        tuningCurve = load('tuningCurveAllOffsets').tuningCurve; %mean rate of spikes per ms.
        SEM = load('sem_values_Tuning').sem_values_Tuning;
        DSI = load(sprintf('Direction-Selectivity-Index-%s',NP.recordingName)).DSI;
        OSI = load(sprintf('Orientation-Tuning-Index-%s',NP.recordingName)).L;


        % Plot the tuning curve with error bars
        fig = figure;
        angles_deg = rad2deg(uDir); % Angles in degrees
        %Convert to spikes per second:
        bar(angles_deg, tuningCurve(u,:).*1000, 'FaceColor', 'k','FaceAlpha',0.5); % Bar plot for spike rates
        hold on;
        errorbar(angles_deg, tuningCurve(u,:).*1000, SEM(u,:).*1000, 'k', 'LineStyle', 'none', 'LineWidth', 1.5); % Error bars
        ylabel('Spikes/sec');
        fig.Position = [  997   266   149   135];

        title(sprintf('OSI = %.2f / DSI = %.2f',OSI(u),DSI(u)))

        if savePlot
            cd(saveDir)
            print(fig, sprintf('%s-MovBall-tuningCurve-U%d.pdf',NP.recordingName,u),'-dpdf', '-r300', '-vector');
        end


    end

ur = ur+1;

end %%%%End for loop over good neurons



end