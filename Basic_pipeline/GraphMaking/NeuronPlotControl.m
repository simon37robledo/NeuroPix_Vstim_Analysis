%%%%% Neuron plot control

% stim = "MB";
% eNeuron = 5;
%%
function NeuronPlotMovingBal(NP,eNeuron,varagin) 


% Parse optional inputs
p = inputParser;

% Add required inputs
addRequired(p, 'NP');
addRequired(p, 'eNeuron');

% Loop through varargin as Name-Value pairs
for i = 1:2:length(varargin)
    name = varargin{i};   % Extract name
    value = varargin{i+1}; % Extract value
    fprintf("x%s = %s\n", name, mat2str(value)); % Print output
end
%% %%%%%%%%%%%%%%%%% load stimulus properties %%%%%%%%%%%%%%%%%%%%%%%%%%%



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
ballFiles = filenames(contains(filenames,"linearlyMovingBall"));


if isempty(ballFiles)
    %disp()
    w= sprintf('No moving ball files where found in %s. Skipping into next experiment.',NP.recordingName);
    warning(w)

end

directions = [];
offsets = [];
sizes = [];
speeds = [];
orientations = [];

j =1;

if size(ballFiles) ~= [0 0]

    for i = ballFiles %Extract stim properties
        ball= load(stimDir+"\"+string(i));

        if ~isempty(find(strcmp(ball.VSMetaData.allPropName,'orientations'))) %Check if orientations are present (grid inside moving object).
            orientations = [orientations cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'orientations'))))];
        else %No orientations property exist
            orientations = [orientations zeros(1,cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials')))))];
        end
        if isempty(orientations) %orientations property exists but is empty
            orientations = [orientations zeros(1,cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials')))))];
        end

        directions = [directions cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'directions'))))];

        offsets = [offsets cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'offsets'))))];

        sizes = [sizes cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballSizes'))))];

        speeds = [speeds cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'speeds'))))];

        interStimStats = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'interTrialDelay'))))*1000;

        j = j+1;
    end
    disp('Visual stats extracted!')
else
    disp('Directory does not exist!');
end

if ~isempty(find(strcmp(ball.VSMetaData.allPropName,'orientationFreq')))
    Freq = [ cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'orientationFreq'))))];
end


%% %%%%%%%%%%%%%%%%%%%% load neuron properties %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = NP.convertPhySorting2tIc(NP.recordingDir);
label = string(p.label');
goodU = p.ic(:,label == 'good');
cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');
GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");


% Load Triggers (diode)
Ordered_stims= strsplit(data.VS_ordered{ex},',');
%containsMB = cellfun(@(x) contains(x,'MB'),Ordered_stims);
containsMB = strcmp(Ordered_stims, 'MB');
ttlInd = find(containsMB);


[stimOn stimOff onSync offSync] = NPdiodeExtract(NP,newDiode,1,"MB",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
[stimOn stimOff onSync offSync] = NPdiodeExtract(NP,0,1,"MB",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A



%When dealing with different speeds, save different stim durs, and create M for each speed
A = [stimOn directions' offsets' sizes' speeds' orientations'];
[C indexS] = sortrows(A,[2 3 4 5 6]);

B = [stimOff directions' offsets' sizes' speeds' orientations'];
[Coff indexSo] = sortrows(B,[2 3 4 5 6]);

stimInter= mean(stimOn(2:end)-stimOff(1:end-1));

includeOnespeed=1;
if includeOnespeed

    A =  A(A(:,5) == max(A(:,5)),:);

    B =  B(B(:,5) == max(B(:,5)),:);

    [C indexS] = sortrows(A,[2 3 4 5 6]);

    [Coff indexSo] = sortrows(B,[2 3 4 5 6]);

    stimOn = A(:,1);
    stimOff = B(:,1);

    %digon  =onDigital(A(:,5) == max(A(:,5)));

    stimInter = interStimStats;

end

stimDur = mean(-stimOn+stimOff);

directimesSorted = C(:,1)';
directimesSortedOff = Coff(:,1)';

try
    NeuronVals = load(sprintf('OneSpeed-%d-NeuronRespCat-%s',max(speeds),NP.recordingName)).NeuronVals;
catch
    NeuronVals = load(sprintf('NeuronRespCat-%s',NP.recordingName)).NeuronVals;
end

uDir = unique(directions);
uSpeed = unique(speeds);
offsetN = length(unique(offsets));
direcN = length(unique(directions));
speedN = length(unique(speeds));
sizeN = length(unique(sizes));
uSize = unique(sizes);
orientN = length(unique(orientations));
nT = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials'))));
trialDivision = nT/(offsetN*direcN*speedN*sizeN*orientN); %Number of trials per unique conditions

orderS = [2 3 4 5;4 2 3 5;5 2 3 4;3 2 4 5];
orderNames = {'dir_off_sizes_speeds';'sizes_dir_off_speeds';'speeds_dir_off_sizes';'off_dir_sizes_speeds'};


posX = squeeze(NeuronVals(:,:,3));
posY = squeeze(NeuronVals(:,:,2));
bin = 1;
bin2 =20;
trialsPerAngle = trialDivision*offsetN*speedN*sizeN*orientN;


directimesSorted = C(:,1)';

preBase = round(stimInter/4);

[Mr] = BuildBurstMatrix(goodU(:,eNeuron),round(p.t/bin2),round((directimesSorted-preBase)/bin2),round((stimDur+preBase*2)/bin2));

[nT,nN,nB] = size(Mr);

Mr2 = [];

pvals = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName),'pvalsResponse').pvalsResponse;

pvals = [eNeuron;pvals(eNeuron)];


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% plot raster MB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sizeN=1;


for u = length(eNeuron)

    j=1;

    if sizeN >1

        mergeTrials = trialDivision;

    else
        mergeTrials = 1;
    end

    if mergeTrials ~= 1 %Merge trials

        for i = 1:mergeTrials:nT

            meanb = mean(squeeze(Mr(i:min(i+mergeTrials-1, end),u,:)),1);

            Mr2(j,:) = meanb;

            j = j+1;

        end
    else
        Mr2 = Mr(:,u,:);
    end

    [nT,nB] = size(Mr2);



    fig =  figure;
    subplot(10,1,[1 6]);
    imagesc(squeeze(Mr2).*(1000/bin2));colormap(flipud(gray(64)));
    %Plot stim start:
    xline(preBase/bin2,'k', LineWidth=1.5)
    %Plot stim end:
    xline(stimDur/bin2+preBase/bin2,'k',LineWidth=1.5)
    ylabel(sprintf('%d Trials',nT));
    title(sprintf('U.%d-Unit-phy-%d',eNeuron,GoodU_or(eNeuron)));

    caxis([0 1])
    dirStart = C(1,2);
    offStart = C(1,3);
    for t = 1:nT*mergeTrials
        if dirStart ~= C(t,2)
            yline(t/mergeTrials+0.5,'k',LineWidth=3);
            dirStart = C(t,2);
        end
        if offStart ~= C(t,3)
            yline(t/mergeTrials+0.5,'k',LineWidth=1);
            offStart = C(t,3);
        end

    end

    hold on


    xticklabels([])
    xlim([0 round(stimDur+preBase*2)/bin2])
    xticks([0 preBase/bin2:300/bin2:(stimDur+preBase*2)/bin2 (round((stimDur+preBase*2)/100)*100)/bin2])
    xticklabels([]);



    %%%%Select 10 random trials where Z-score is above the average

    SelectRand = 0;

    if SelectRand

        ZscoreRaster = load(sprintf('ZscoreRaster-%d-%s',N_bootstrap,NP.recordingName)).ZscoreRaster;

        mZscoreRaster = mean(ZscoreRaster(:,eNeuron,:),3);

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


        [maxResp,maxRespIn]= max(NeuronVals(eNeuron,:,1));

        maxRespIn = maxRespIn-1;

        trials = maxRespIn*trialDivision+1:maxRespIn*trialDivision + trialDivision;

        colors = turbo(trialDivision);

        y1 = maxRespIn*trialDivision + trialDivision+1;
        y2 = maxRespIn*trialDivision+2;

        patch([preBase/bin2 (stimDur+preBase)/bin2 (stimDur+preBase)/bin2 preBase/bin2],...
            [y2 y2 y1 y1],...
            'b','FaceAlpha',0.3,'EdgeColor','none')

    end

    %% %%%% Plot PSTH
    %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    subplot(10,1,[7 8])

    MRhist = BuildBurstMatrix(goodU(:,eNeuron),round(p.t),round((directimesSorted(trials)-preBase)),round((stimDur+preBase*2)));

    [nT2,nN2,nB2]= size(MRhist);

    MRhist = squeeze(MRhist);

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
    b.FaceColor = 'k';
    b.FaceAlpha = 0.5;
    b.MarkerEdgeColor = "none";
    ylabel('Firing Rate (Hz)');
    xlim([0 round((stimDur+preBase*2)/100)*100])

    xticks([0 preBase:300:(stimDur+preBase*2) round((stimDur+preBase*2)/100)*100])
    xticklabels([])
    xline([preBase stimDur+preBase],'LineWidth',1.5)


    %set(gca, 'XLim', [-0.5, 1]); % Match raster plot x-axis


    %%%%PLot raw data several trials one
    %%%%channel%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    chan = goodU(1,eNeuron);

    startTimes = directimesSorted(trials)-preBase;

    window = round(stimDur+preBase*2);

    freq = "AP"; %or "LFP"

    type = "line"; %or heatmap

    subplot(10,1,[9 10])

    PlotRawData(fig,NP,chan,startTimes,window,freq,type,preBase);

    xticklabels([]) %%%plot example PV139-local neurons
    xlabel([])
    %xlim([0 round(stimDur+preBase*2)])

    xticklabels([-(preBase) 0:300:round((stimDur/100))*100 round((stimDur/100))*100 + preBase]./1000)
    fig.Position = [680    42   560   954];
    cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
    print(gcf, sprintf('%s-MovBall-SelectedTrials-eNeuron-%d.pdf',NP.recordingName,eNeuron), '-dpdf', '-r300', '-vector');


    %% %%%%%%%%%Plot receptive field per direction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cd(NP.recordingDir)

    %Parameters
    %check receptive field neurons first
    respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;

    respU = find(respNeuronsMB<0.05);

    ru = find(respU == eNeuron);
    %
    reduceFactor =20;
    eye_to_monitor_distance = 21.5; % Distance from eye to monitor in cm
    pixel_size = 33/(1080/reduceFactor); % Size of one pixel in cm (e.g., 25 micrometers)
    monitor_resolution = [redCoorX, redCoorY]; % Width and height in pixels
    [theta_x,theta_y] = pixels2eyeDegrees(eye_to_monitor_distance,pixel_size,monitor_resolution);

    coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));
    redCoorX = round(coorRect(3)/reduceFactor);
    redCoorY = round(coorRect(4)/reduceFactor);

    theta_x = theta_x(:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY);



    if noEyeMoves
        q=1;
        RF{1} = load(sprintf('NEM-RFuSTDirSizeFilt-Q%d-%s',q,NP.recordingName)).RFuSTDirSizeFilt; %Size and dir
        q=2;
        RF{2} = load(sprintf('NEM-RFuSTDirSizeFilt-Q%d-%s',q,NP.recordingName)).RFuSTDirSizeFilt; %Size and dir
    else
        q=1;
        RFuSTDirSizeFilt = load(sprintf('RFuSTDirSizeFilt-%s',NP.recordingName)).RFuSTDirSizeFilt; %Size and dir
    end

    if size(RFuSTDirSizeFilt,2) >1 %If there is more than on size %%Solve for 2 quadrants
        %%Select size closest to 120
        [minS indx] = min(abs(unique(sizeV)-120));
        RFuSTDirFilt = squeeze(RFuSTDirSizeFilt(:,indx,:,:,:));
    else
        RFuSTDirFilt = squeeze(RFuSTDirSizeFilt);
    end
    %
    for q = numel(RF)

        RFuSTDirFilt = squeeze(RF{q});

        figure;
        DirLayout=tiledlayout(direcN/2,direcN/2,"TileSpacing","tight","Padding","none");


        for d = 1:direcN

            nexttile(DirLayout);

            if d ==1 || d==3
                imagesc(rot90(squeeze(RFuSTDirFilt(d,:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY,ru)),2));caxis([0 max(RFuSTDirFilt(:,:,:,ru),[],'all')]);c = colorbar;
            else
                imagesc((squeeze(RFuSTDirFilt(d,:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY,ru))));caxis([0 max(RFuSTDirFilt(:,:,:,ru),[],'all')]);c = colorbar;
            end
            colormap('turbo')
            title(string(uDir(d)))
            title(c,'Z-score')
            %xlim([(redCoorX-redCoorY)/2 (redCoorX-redCoorY)/2+redCoorY])
            xt = xticks;%(linspace((redCoorX-redCoorY)/2,(redCoorX-redCoorY)/2+redCoorY,offsetN*2));
            xticklabels(round(theta_x(1,xt)))
            yt = yticks;
            yticklabels(round(theta_y(yt,1)))
            xlabel('X degrees')
            ylabel('Y degrees')

            axis equal tight


        end

        title(DirLayout,sprintf('%d-Q-%d',respU(ru),q),'FontSize',10)

        %     cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
        %     print(gcf, sprintf('%s-MovBall-ReceptiveField-q-%deNeuron-%d.pdf',NP.recordingName,q,eNeuron), '-dpdf', '-r300', '-vector');

        %%%%%%%%%

        spikeSum = spikeSums{q};

        eNspkS = squeeze(spikeSum(:,ru,:));
        figure;imagesc(eNspkS);colormap(flipud(gray(64)));
        %     rowsWithNaN = find(any(isnan(eNspkS), 2));
        %     yline(rowsWithNaN,'g')
        yline(trialDivision*sizeN:trialDivision*sizeN:size(spikeSums{q},1));
        yline(trialDivision*offsetN*sizeN:trialDivision*offsetN*sizeN:size(spikeSums{q},1)-1,'r','LineWidth',5)
        title(string(q))

    end %%%%%%%%%%%%%%%%


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%Plot eye movements colorcoded for these selected trials:

    %%%%%% Construct stimType matrix for eye movement plot.
    stimType = zeros(length(A),5); %3 plus number of example neuron
    stimType(:,1) = A(:,1);
    stimType(:,2) = A(:,1)+stimDur;
    stimType(:,3) = A(:,2);
    stimType(:,4) = A(:,3);

    %Get response strenght of specific neurons and save it in stimType
    [MrNoSort] = BuildBurstMatrix(goodU(:,u),round(p.t/bin),round((stimOn'-preBase)/bin),round((stimDur+preBase*2)/bin)); %response matrix
    ResponseStrength= mean(MrNoSort(:,u,round(preBase/bin):round((preBase+stimDur)/bin)),3); %For PV35_3
    stimType(:,end) = ResponseStrength;

    timesnips = [directimesSorted(trials);directimesSorted(trials)+stimDur];

    %EyePositionAnalysis(NP,data.Eye_video_dir{ex},21,1,0,0,timesnips,1,0)

    EyePositionAnalysis(NP,data.Eye_video_dir{ex},21,1,0,0,0,1,1)

    cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\1stFigure')
    print(gcf, sprintf('%s-MovBall-EyeMovs-SelectedTrials-eNeuron-%d.pdf',NP.recordingName,eNeuron), '-dpdf', '-r300', '-vector');

    %%%Plot ellipses in each color corresponding to the mean point of
    %%%each trial.

    %%Mark trials and plot raw data of trials



    %
    %         close

end


end