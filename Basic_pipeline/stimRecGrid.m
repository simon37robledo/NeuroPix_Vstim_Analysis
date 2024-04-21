%% PROCESSING RECTANGLE GRID STIMULUS

%% Phy commmand: phy template-gui params.py
%Order of batch= FFF, movBall, rectGrid, RectNoise, Gratings.
%% Phy commmand: phy template-gui params.py
%Order of batch= FFF, movBall, rectGrid, RectNoise, Gratings.
ins = 1:2;
basic_path = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV102';
exp = 'PV102_experiment_18_7_23';

ttl_index =3;
in_ind =1;

zScoreRate= cell(1,length(ins));

allDepth =  cell(1,length(ins));

allXdist =  cell(1,length(ins));

spkRateMinus = cell(1,length(ins));

spkRateRmean = cell(1,length(ins));
spkRateBmean = cell(1,length(ins));
spkRateAmean = cell(1,length(ins));


animal = "PV102";
rasters = 1;
rastersFig = 1;
heatmaps = 0;
heatmap_fig = 0;
rawTrace = 0;
rawTraceFig = 0;
NumBatch = 1; 
spikes = 1;
rangeMV = [-100 100];

%SA6
% InsAngles = [88 88 88];
% InsDepths = [3931 3955 3992];

% %PV102
 InsAngles = [88 88];
 InsDepths = [3909 3931];

% %PV67
%InsAngles = [88 88 88 88];
%InsDepths = [3956 3907 4000 3934];

%%
for in = 2 %Start iterating across insertions


    path = convertStringsToChars(basic_path+string(filesep)+exp+string(filesep)+"Insertion"+in+string(filesep)+"catgt_"+exp+"_"+in+"_g0");

    NP = NPAPRecording(path);
  

    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

    cluster_group = readtable(string(NP.recordingDir) + "\cluster_group.tsv",  "FileType","text",'Delimiter', '\t');

    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");


    AngleInser = InsAngles(in_ind);
    InserDepth = InsDepths(in_ind);
    ShankDist = zeros(1, length(GoodU_or)); %distance of unit from surface along the shank
    verticalDepth = zeros(1, length(GoodU_or)); %depth of unit along vertical axis
    XDist = zeros(1, length(GoodU_or)); %X distance of unit from insertion
    Zscore = zeros(1, length(GoodU_or));
    minusSpkrt = zeros(1, length(GoodU_or));
    spkRateA = zeros(1, length(GoodU_or));
    spkRateBM = zeros(1, length(GoodU_or));
    spkRateRM = zeros(1, length(GoodU_or));


  

    [Ttrigger,chNumberT]=NP.getTrigger();

    Sstart = [Ttrigger{1}];
    Send = [Ttrigger{2}];
    onset = [Ttrigger{3}];
    offset = [Ttrigger{4}];

    stimOn = [];
    stimOff = [];
    stimInter = [];
    ttlNum = zeros(1,length(ttl_index));
    j =1;

    for i=ttl_index
        stimUp = onset(onset > Sstart(i) & onset < Send(i));
        stimOn = [stimOn stimUp]; %general

        stimDown = offset(offset > Sstart(i) & offset < Send(i));
        stimOff = [stimOff stimDown]; %general

        stimInter = [stimInter mean(stimOff-stimOn)];

        %ttlNum(j) = length(stim(1:2:end)); %sanity check to see how many stimulus presentations there are per round

        j = j+1;
    end

    stimDur = stimOff(1) - stimOn(1);

    stimInter = mean(stimOn(2:end)-stimOff(1:end-1));

    stimOn = stimOn'; %transform into horizontal vector.



    %--------------------------------------------------------------------------------------------------------------------------------%
    % Matlab visual stimulus statistics
    stimDir =convertStringsToChars(basic_path+"\"+exp+"\Insertion"+in);

    file = dir (stimDir);
    filenames = {file.name};
    recFiles = filenames(contains(filenames,"rectGrid"));
    seqMatrix = [];

    if size(recFiles) ~= [0 0]

        j =1;
        for i = recFiles
            rec= load(stimDir+"\"+string(i));


            seqMatrix = [seqMatrix cell2mat(rec.VSMetaData.allPropVal(16))];

            positionsMatrix = [cell2mat(rec.VSMetaData.allPropVal(14)),cell2mat(rec.VSMetaData.allPropVal(15))];

            stimDurStats = cell2mat(rec.VSMetaData.allPropVal(38))*1000;
            interStimStats = cell2mat(rec.VSMetaData.allPropVal(28))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    %--------------------------------------------------------------------------------------------------------------------------------%
    % Plot raw data according to rectangle positions that is presented
    if rawTrace ==1

        pathFig = '\\132.66.45.127\data\Large_scale_mapping_NP\Immobilized_exp\PV102\PV102_experiment_18_7_23\Insertion1\Figures\RectGrid';

        cd(pathFig)


        u = 1:16:384;

        FH=filterData(NP.samplingFrequency);
        FH.highPassPassCutoff=100;
        FH.highPassStopCutoff=80;
        FH.lowPassPassCutoff=1800;
        FH.lowPassStopCutoff=2000;
        FH.attenuationInLowpass=20;
        FH.attenuationInHighpass=20;
        FH=FH.designBandPass;

        window = stimDurStats+interStimStats;

        HP = 1;



        for u = u


            %Raw traces:

            channel = u;


            trace = zeros(numel(unique(seqMatrix)),cell2mat(rec.VSMetaData.allPropVal(29))*NumBatch,round(window*NP.samplingFrequency/1000));

            for i = unique(seqMatrix)
                start_times = stimOn(seqMatrix == i);

                %Option 1 = mean of x channels
                %d = NP.getData(channel:(channel+15),start_times-50,stimDur);OR
                %get mean of 15 channels
                %trace(i,:,:) = mean(d, 1); %Mean across channels

                %Option 2 = select one channel every x channels and do HP filter or not
                d = NP.getData(channel,start_times-stimDurStats/2,window);
                if HP == 0
                    trace(i,:,:) = d; %Select one channel
                    type = 'LFP';
                else
                    trace(i,:,:) = FH.getFilteredData(d); %filter data
                    type = 'HP';
                end

            end


            En =  fliplr(reshape(sort(unique(seqMatrix),"descend"),[sqrt(max(seqMatrix)), sqrt(max(seqMatrix))])');

            title_stim = sprintf('Insertion-%d-RectGrid-Raw-channel-%d-%s', in, channel, type);

            [nChan,ntrials,nTimes]=size(d);

            h = figure('units','normalized','outerposition',[0 0 1 1]);
            hplot = activityMultiTracePhysicalSpacePlot([], unique(seqMatrix),trace,En,'scaling','std');
            title(sprintf('Raw traces: Channel %d',channel))
            lim = get(gca, 'XLim'); %Get strange limits of divided plot
            xline((1.5 : 2*(sqrt(max(seqMatrix)))+0.5)*lim(1), '--k');
            xticks((1.5 : 2*(sqrt(max(seqMatrix)))+0.5)*lim(1));
            xticklabels(repmat([0 window/2],1,sqrt(max(seqMatrix))));
            xlabel('[Ms] from stimulus onset');
            yticks((1 : 2*(sqrt(max(seqMatrix)))+0.05)*lim(1));
            yticklabels(repmat([ntrials ntrials/2],1,sqrt(max(seqMatrix))));
            ylabel(sprintf('%d trials per x-y location', ntrials));
            print(sprintf('%s.png',title_stim),'-dpng');

            close(h)

            clear h hplot trace


        end

    end
    %--------------------------------------------------------------------------------------------------------------------------------%
    % Load spike sorting and Phy results

    % Convert into format to generate raster plots:
    p = NP.convertPhySorting2tIc(NP.recordingDir);

    %Select good units
    label = string(p.label');

    goodU = p.ic(:,label == 'good');
    amp = p.neuronAmp(label=='good');

    disp("Phy extracted!")

   
    %--------------------------------------------------------------------------------------------------------------------------------%
    % Build raster plots per unit
    bin=20;
    win=stimDurStats+interStimStats;
    start = interStimStats/2;

    %image saving path:
    pathF = '\\132.66.45.127\data\Large_scale_mapping_NP\Immobilized_exp\PV102\PV102_experiment_18_7_23\Insertion1\Figures\RectGrid';
    cd(pathF)


    ps = 200; %peristimulus time;

    TileNum = length(unique(seqMatrix));

    allM = zeros(length(GoodU_or),sqrt(TileNum),sqrt(TileNum));


    for u = 1:length(GoodU_or)

        u = 184;

        %Position per unit:

        cluster_info.ch(cluster_info.ch ==0) = 1; %Change channel 0 to 1.

        ch = cluster_info.ch(cluster_info.cluster_id == GoodU_or(u));

        title_stim = sprintf('Insertion-%d-FFF-Unit-%d-channel-#%d', in, GoodU_or(u), cluster_info.ch(cluster_info.cluster_id == GoodU_or(u)));


        ShankDist(u) = InserDepth - (NP.chLayoutPositions(2,ch));

        verticalDepth(u) = sin(deg2rad(AngleInser))*ShankDist(u); %depth of unit along vertical axis

        XDist(u) = cos(deg2rad(AngleInser))*ShankDist(u); %X distance of unit from insertion

        %Title:
        title_stim = sprintf('Insertion-%d-RectGrid-Unit-%d-channel-#%d', in, GoodU_or(u), cluster_info.ch(cluster_info.cluster_id == GoodU_or(u)));


        t = tiledlayout(sqrt(max(seqMatrix)), sqrt(max(seqMatrix)),'TileSpacing','tight');

        HM = reshape( 1:length(unique(seqMatrix)), 10, 10) .';
        HZ = zeros(10);
        HF =  zeros(10);


        if rasters == 1
            for i = 1:max(seqMatrix)
                %Build raster
                bin = 60;
                start_times = stimOn(seqMatrix == i)-start;
                [M]=BuildBurstMatrix(goodU,round(p.t/bin),round(start_times/bin),round(win/bin));
                [nTrials,nNeurons,nTimes]=size(M);

                if rastersFig ==1

                    nexttile
                    imagesc((1:nTimes)*bin,1:nTrials,squeeze(M(:,u,:)));colormap(flipud(gray(64)));
                    xline(start, '-', LineWidth=3, Color="#77AC30");
                    xline(stimDurStats+start, '-', LineWidth=3, Color="#0072BD");
                    xticks([start stimDurStats+start win]);
                    caxis([0 1]);
                    set(gca,'YTick',[]);
                    set(gca,'YTickLabel',[]);


                    if i < max(seqMatrix)-sqrt(max(seqMatrix)-1)
                        set(gca,'XTickLabel',[]);

                    end

                end

                start_times = stimOn(seqMatrix == i);

                %Rate over miliseconds:
                bin=stimDur+ps;

                [Mps] = BuildBurstMatrix(goodU,round(p.t/bin),round(start_times/bin),round((stimDur+ps)/bin)); %response

                AllSpks = squeeze(M(:,u,:));

                ResSpks = Mps(:,u);

                HZ(HM == seqMatrix(i)) = mean(ResSpks, 'all');

            end
            fig = gcf;
            set(fig, 'Color', 'w');
            colorbar

    % Set the color of the figure and axes to black
            title(t,title_stim)
            ylabel(t,sprintf('%d trials',nTrials))
        end

        if heatmaps ==1
        bin = stimInter/2;
        [Mbs]=BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn-stimInter/2)/bin),(stimInter/2)/bin); %baseline in rate

        BaseSpks = Mbs(:,u);

        HBm = mean(BaseSpks, 'all');

        HBstd = std(BaseSpks,0, 'all');

        Hn = (HZ-HBm)/HBstd; %Normalization: Z score

        Hn(isnan(Hn)|isinf(Hn))=0;

        LB=flipud(lbmap(256,'BrownBlue'));

        Hrez = imresize(Hn, [7 7], 'bilinear');
        figure
        imagesc(Hrez);
        %plot it
        title(title_stim);
        colormap(LB);
        caxis([-2 2]);
        % caxis([-max(abs(Hn(:)))-0.1 max(abs(Hn(:)))+0.1]);
        ax = gca;
        %set(gca,'color',0*[1 1 1]);
        set(gca,'XTick',[],'YTick',[]);
        axis equal
        axis tight
        set(gcf,'PaperPositionMode','auto')
        cb = colorbar();
        cb.Ticks = [-2, 1, 0, 1, 2];
        print(sprintf('%s-MovBall-Positions-mean',title_stim),'-dpdf','-vector');

        fig = gcf;
        set(fig, 'Color', 'w');

        allM(u,:,:) = Hn;
        close all
        end

        %%% Calculate Z-score:

        %Baseline
        stim = stimOn;

        [Mb]=BuildBurstMatrix(goodU,round(p.t),round(stim-stimInter),round(stimInter));

        spkRateB = sum(squeeze(Mb(:,u,:)),2)/round(stimInter/1000);

        spkRateBM(u) = mean(spkRateB);

        %Respones

        [Mr]=BuildBurstMatrix(goodU,round(p.t),round(stim),round(stimDur));


        spkRateR = sum(squeeze(Mr(:,u,:)),2)/round(stimDur/1000);

        spkRateRM(u) = mean(spkRateR);

        %All
        stimRoundDur = stimOff(end)+ stimInter- stim(1);

        [Ma] = BuildBurstMatrix(goodU,round(p.t),round(stim(1)),round(stimRoundDur));
        
        spkRateA(u) = sum(squeeze(Ma(:,u,:)),'all')/round(stimRoundDur/1000);

        Zscore(u) = (spkRateRM(u)-spkRateBM(u))/std(spkRateB');

        minusSpkrt(u) = mean(spkRateR-spkRateB);

    end
        allDepth{in_ind} = {verticalDepth};

        allXdist{in_ind} = {XDist};

        zScoreRate{in_ind} = {Zscore};

        spkRateMinus{in_ind} = {minusSpkrt};

        spkRateRmean{in_ind} = {spkRateRM};

        spkRateBmean{in_ind} = {spkRateBM};

        spkRateAmean{in_ind} = {spkRateA};

        in_ind = in_ind+1;

end

save(sprintf('%s-RG-spkR_zscore',exp),"zScoreRate")
save(sprintf('%s-unit-depths',exp),"allDepth")
save(sprintf('%s-unit-x-dist',exp),"allXdist")
save(sprintf('%s-RG-spkR_minus',exp),"spkRateMinus")

save(sprintf('%s-RG-spkRateRmean',exp),"spkRateRmean")
save(sprintf('%s-RG-spkRateBmean',exp),"spkRateBmean")
save(sprintf('%s-RG-spkRateAmean',exp),"spkRateAmean")


%%

allMm = (squeeze(mean(allM,1,"omitnan")));

imagesc(allMm); %plot it
colormap(LB);
caxis([-5.63193668669677 5.63193668669677]);
ax = gca;
set(gca,'XTick',[],'YTick',[]);
axis equal
axis tight
colorbar()
%ax.CLim = [0  maxVal];
%colorbar('gray');
title(sprintf('%s-rectGrid-mean',title_stim))
%     xlabel(t,'X')
%     ylabel(t,'Y');

%set(gcf,'PaperPositionMode','auto')
print(sprintf('%s-rectGrid-mean',title_stim),'-dpdf','-vector');



