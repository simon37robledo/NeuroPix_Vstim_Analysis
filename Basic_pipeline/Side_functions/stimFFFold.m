%% PROCESSING FULL FIELD FLASH STIMULUS


%% Phy commmand: phy template-gui params.py
%Order of batch= FFF, movBall, rectGrid, RectNoise, Gratings.
ins = 1:2;
basic_path = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV102';
exp = 'PV102_experiment_18_7_23';
animal = "PV102";
in_ind =1; %initialize insertion index. 
ttl_index =1;
in_ind =1;

zScoreRateB= cell(1,length(ins));

zScoreRateA= cell(1,length(ins));

allDepth =  cell(1,length(ins));

allXdist =  cell(1,length(ins));

spkRateMinus = cell(1,length(ins));

spkRateRmean = cell(1,length(ins));
spkRateBmean = cell(1,length(ins));
spkRateAmean = cell(1,length(ins));
spkRateAtmean = cell(1,length(ins));


rasters = 1;
rastersFig = 1;
heatmaps = 0;
heatmap_fig = 0;
rawTrace = 0;
rawTraceFig = 0;
NumBatch = 1; 
spikes = 1;
rangeMV = [-200 200];

%SA6
% InsAngles = [88 88 88];
% InsDepths = [3931 3955 3992];

%PV102
InsAngles = [88 88];
InsDepths = [3909 3931];

% % %PV67
% InsAngles = [88 88 88 88];
% InsDepths = [3956 3907 4000 3934];


%% Iterate through insertions

for in = 2


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
    ZscoreB = zeros(1, length(GoodU_or));
    minusSpkrt = zeros(1, length(GoodU_or));
    spkRateA = zeros(1, length(GoodU_or));
    spkRateBM = zeros(1, length(GoodU_or));
    spkRateRM = zeros(1, length(GoodU_or));
    ZscoreA = zeros(1, length(GoodU_or));
    spkRateAt = zeros(1, length(GoodU_or));


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

 stimInter =  mean(stimOn(2:end)-stimOff(1:end-1));

 stimOn = stimOn'; %transform into horizontal vector. 


%--------------------------------------------------------------------------------------------------------------------------------%
% Raw traces channel examples

pathF = convertStringsToChars(basic_path+"\"+exp+"\Insertion"+in+"\Figures\FFF");

cd(pathF)

if rawTrace == 1
ch = 1:16:384;
window = stimDur+stimInter;
stimStart = stimInter/2;

d = NP.getData(ch,stimOn-stimStart,window);

FH=filterData(NP.samplingFrequency);
FH.highPassPassCutoff=100;
FH.highPassStopCutoff=80;
FH.lowPassPassCutoff=1800;
FH.lowPassStopCutoff=2000;
FH.attenuationInLowpass=20;
FH.attenuationInHighpass=20;
FH=FH.designBandPass;

for i = 1:length(ch)
    h = figure('units','normalized','outerposition',[0 0 1 1]);
    
    if HPfilter == 1
    di = FH.getFilteredData(d(i,:,:)); f=10; %for HP filter

    else
    di = d(i,:,:); f=5;
    end

    di = squeeze(di);
    s = size(di);
    
    %Plot raw trace.
    if rawTraceGig ==1
    title_stim = sprintf('Insertion-%d-FFF-Raw-channel-%d', in, ch(i));
    offset = f*std(di); % Make this big enough to prevent overlap
    offset_vector = (offset:offset:s(1)*offset)';
    M_plus_offset = bsxfun(@plus,di,offset_vector);
    plot(M_plus_offset')
    xline(stimStart*NP.samplingFrequency/1000);
    xline((stimStart+stimDur)*NP.samplingFrequency/1000);
    xticks([0:100:window]*NP.samplingFrequency/1000);
    xticklabels(0:100:window);
    set(gca,'ytick',[]);
    title(sprintf('FFF Raw traces: Channel %d',ch(i)))
    xlabel('[Ms]');
    ylabel(sprintf('%d Trials', length(stimOn)));
    print(sprintf('%s.png',title_stim),'-dpng');  
    close(h)
    clear h 
    end

end
end

%--------------------------------------------------------------------------------------------------------------------------------%
%LFP as heatmap:

if heatmaps ==1

    chanNUM = numel(NP.channelNumbers(1:end-1));
    
   channels = 197;

    lfd = NP.getData(channels,stimOn-stimInter/2,stimDur+stimInter);

    lfd = bsxfun(@minus,lfd,mean(lfd,3)); % Substract the mean from each sample.
    %lfd = lfd(1:end ~= 54,:,:);-
    [nCH,trials,samples] = size(lfd);

    if length(channels)~=1
        y = 1:nCH;
        imageLF = squeeze(mean(lfd,2));
        ylab = 'Channels';
    else
        y = 1:nTrials;
        imageLF = squeeze(lfd);
        ylab = sprintf('%d Trials',length(stimOn));
        yt=1:length(stimOn);
    end

if heatmap_fig ==1
    title_stim = (sprintf('In. %d FFF channels %d %d',in, max(channels),min(channels)));
    h = figure;
    imagesc(1:samples/round(NP.samplingFrequency),y,imageLF);colormap(parula(64));ylabel(ylab);xlabel('Time [ms]');
    title(title_stim);
    yticks([])
    cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Mean LFP (mV)');
    line([stimInter/2 stimInter/2],ylim,'color','b','LineWidth',4);
    line([stimInter/2+stimDur stimInter/2+stimDur],ylim,'color','r','LineWidth',4);
    set(gca,'YDir','normal'); %reverse image

    set(h, 'Color', 'w');
    saveas(h,sprintf('%s.png',title_stim));
    clear h
end

end
%--------------------------------------------------------------------------------------------------------------------------------%

% Load spike sorting and Phy results

if spikes == 1
    % Convert into format to generate raster plots:
    p = NP.convertPhySorting2tIc(NP.recordingDir);

    %Select good units
    label = string(p.label');

    goodU = p.ic(:,label == 'good');

    disp("Phy extracted!")

    % %--------------------------------------------------------------------------------------------------------------------------------%
    %BUILD RASTER PLOT
    bin=20;
    win=stimDur+stimInter*2;
    start = stimInter;


    tDiode = "False"; %use diode times?

    if tDiode == "True"
        stim = tDiodeOr;
    else
        stim = stimOn;
    end

    e =1; %error index

    errorU = [];

    for u = 1:length(GoodU_or)

        u = 184;
        %Position per unit:

        cluster_info.ch(cluster_info.ch ==0) = 1; %Change channel 0 to 1.

        ch = cluster_info.ch(cluster_info.cluster_id == GoodU_or(u));


        ShankDist(u) = InserDepth - (NP.chLayoutPositions(2,ch));

        verticalDepth(u) = sin(deg2rad(AngleInser))*ShankDist(u); %depth of unit along vertical axis

        XDist(u) = cos(deg2rad(AngleInser))*ShankDist(u); %X distance of unit from insertion

        %Title:
        title_stim = sprintf('Insertion-%d-Directions-Unit(%d)-%d-channel-#%d', in, GoodU_or(u),u, cluster_info.ch(cluster_info.cluster_id == GoodU_or(u)));

        if rasters ==1
            %Build burst matrix.
            stimS = stim- start;

            [M]=BuildBurstMatrix(goodU,round(p.t/bin),round(stimS/bin),round(win/bin));

            [nTrials,nNeurons,nTimes]=size(M);

            if rastersFig == 1

                %Plot average accros trials where Units are the y axis.
                h = figure();
                imagesc((1:nTimes)*bin,1:nTrials,squeeze(M(:,u,:)));colormap(flipud(gray(64)));ylabel('Trials');xlabel('Time (ms)');
                title(sprintf('%s -%s trials.',title_stim, string(nTrials)));
                xline(start, '-g', 'Start','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
                xline((start+stimDur), '-b', 'End', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
                %yticklabels(1:nTrials);
                try
                    print(h,title_stim, '-dpng');
                catch ME
                    fprintf('Error occurred when saving unit %d:  %s\n', GoodU_or(u), ME.message);
                    errorU = [errorU u];
                    e = e+1;
                    continue; % Skip to the next iteration
                end

                clearvars h
                close all
            end
        end


        %%% Calculate Z-score:

        %Baseline

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

        %All per tiral

        [Mat] = BuildBurstMatrix(goodU,round(p.t),round(stim),round(stimDur+stimInter));
        
        spkRateAtrials = sum(squeeze(Mat(:,u,:)),2)/round((stimDur+stimInter)/1000);

        spkRateAt(u) = mean(spkRateAtrials);

        ZscoreA(u) = (spkRateRM(u)-spkRateAt(u))/std(spkRateAtrials');

        ZscoreB(u) = (spkRateRM(u)-spkRateBM(u))/std(spkRateB');

        minusSpkrt(u) = mean(spkRateR-spkRateB);
    end


    allDepth{in_ind} = {verticalDepth};

    allXdist{in_ind} = {XDist};

    zScoreRateB{in_ind} = {ZscoreB};
    
    zScoreRateA{in_ind} = {ZscoreA};

    spkRateMinus{in_ind} = {minusSpkrt};

    spkRateRmean{in_ind} = {spkRateRM};
    
    spkRateBmean{in_ind} = {spkRateBM};
    
    spkRateAmean{in_ind} = {spkRateA};

    spkRateAtmean{in_ind} = {spkRateAt};

    in_ind = in_ind+1;


end

end
cd(basic_path+string(filesep)+exp)
save(sprintf('%s-NonSaved-units',exp),"errorU")

save(sprintf('%s-FFF-spkR_zscoreB',exp),"zScoreRateB")
save(sprintf('%s-FFF-spkR_zscoreA',exp),"zScoreRateA")

save(sprintf('%s-unit-depths',exp),"allDepth")
save(sprintf('%s-unit-x-dist',exp),"allXdist")

save(sprintf('%s-FFF-spkR_minus',exp),"spkRateMinus")
save(sprintf('%s-FFF-spkRateRmean',exp),"spkRateRmean")
save(sprintf('%s-FFF-spkRateBmean',exp),"spkRateBmean")
save(sprintf('%s-FFF-spkRateAmean',exp),"spkRateAmean")
save(sprintf('%s-FFF-spkRateAmeanT',exp),"spkRateAtmean")




%% Plot insertions

cd(basic_path+string(filesep)+exp)

animalZSA = load(sprintf('%s-FFF-spkR_zscoreA',exp));

animalZSB = load(sprintf('%s-FFF-spkR_zscoreB',exp));

animalZSb = animalZSB.zScoreRateB;

animalZSa = animalZSA.zScoreRateA;

animalM = load(sprintf('%s-FFF-spkR_minus',exp));
animalM = animalM.spkRateMinus;

AnimalR =  load(sprintf('%s-FFF-spkRateRmean.mat',exp));
animalR = AnimalR.spkRateRmean;

AnimalB = load(sprintf('%s-FFF-spkRateBmean.mat',exp)); 
animalB = AnimalB.spkRateBmean;

AnimalAll = load(sprintf('%s-FFF-spkRateAmean.mat',exp));
animalAll =  AnimalAll.spkRateAmean; 

allDepth = load(sprintf('%s-unit-depths',exp));


%%
figure()

InsN = 2;


if InsN > 2
    t = tiledlayout(round(InsN/2),round(InsN/2),'TileSpacing','tight');
else
    t = tiledlayout(InsN,1,'TileSpacing','tight');
end


maxZ = 0;

minZ = inf;

for i = 1:InsN %insertions



    DepthMat = cell2mat(allDepth.allDepth{i});
    zscoreMatB = cell2mat(animalZSb{i});
    zscoreMatA = cell2mat(animalZSa{i});
    minus = cell2mat(animalM{i});
    NormB = minus./cell2mat(animalB{i});
    NormAll = cell2mat(animalR{i})./cell2mat(animalAll{i});

    colors = rand(length(DepthMat),3);
    %colors(colors(:, 1) < 0.2, :) = [1, 0, 0];

    currentZM = max(zscoreMatB(~isinf(zscoreMatB)));

    currentZm = min(zscoreMatB); 

    if currentZM > maxZ
        maxZ = currentZM;
    end
    
    if currentZm < minZ
        minZ = currentZm;
    end
    
    yAxis = round(DepthMat)*-1;

    transp = 0.1;
    tiles(i) = nexttile;
    signif = zscoreMatB > 2 & zscoreMatB < inf;

    scatter(minus(signif),yAxis(signif),[], zscoreMatB(signif) , "filled")
    hold on
    title(sprintf("%s: Insertion %d - FFF",animal,i))
    scatter(minus(~signif),yAxis(~signif), [], zscoreMatB(~signif), "filled", MarkerFaceAlpha=transp)
    hold off
    hcb=colorbar;
    hcb.Title.String = "Z-score";
    colormap(flipud(gray));
    caxis([floor(minZ), ceil(maxZ)]);
    

end

xlabel(t,'SpkR (response - baseline)');
ylabel(t,'Depth (um)');
linkaxes(tiles);
linkprop(findall(gcf, 'Type', 'ColorBar'), 'Limits');

print(sprintf('%s-FFF',exp), '-dpng')

GoodU_or(signif)
