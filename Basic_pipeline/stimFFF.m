%% PROCESSING FULL FIELD FLASH STIMULUS

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%Optionall
RawData = 0;
plotexamplesFF =0;
newTIC = 0;
ex=1;
ZscoresDo=0;
plotin3D=0;
processLFP=0;
Shuffling_baseline = 1;
repeatShuff = 1;

plotRasters = 0;

%%
% Iterate through experiments (insertions and animals) in excel file
for ex =GoodRecordings%size(data,1)
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

     %2. Extract moving ball statistics
    patternIndex = strfind(string(NP.recordingDir), "\catgt");

    endIndex = patternIndex(1)-1;
    stimDir = string(NP.recordingDir);
    stimDir = extractBetween(stimDir,1,endIndex);

    file = dir (stimDir);
    filenames = {file.name};

    file = dir (stimDir);
    filenames = {file.name};
    FFFfiles = filenames(contains(filenames,"fullFieldFlash"));

    if isempty(FFFfiles)
        %disp()
        w= sprintf('No full field flash files where found in %s. Skipping into next experiment.',NP.recordingName);
        warning(w)
        continue
    end

    %3. Load Triggers (diode)
    Ordered_stims= strsplit(data.VS_ordered{ex},',');
    containsFF = cellfun(@(x) contains(x,'FFF'),Ordered_stims);
    ttlInd = find(containsFF);

    [stimOn stimOff] =  NPdiodeExtract(NP,0,0,"FFF",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff] =  NPdiodeExtract(NP,0,0,"FFF",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A

    stimDur = mean(stimOff-stimOn);
    stimInter= mean(stimOn(2:end)-stimOff(1:end-1));
    %%% Raw traces
    %Channel 219 - insertion1 PV67

    p = NP.convertPhySorting2tIc(NP.recordingDir);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');

    bin = 1;
    preBase = round(stimInter/2);
    stimOn = stimOn';

    duration = 300; %Ms where off response can be found

    Mr = BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn)/bin),round((stimDur+duration)/bin));

   
    [nT,nN,nB] = size(Mr);

    responseWindow = 20:floor(stimDur)+duration;

%     artN = zeros(nT,1,nB);
% 
%     artN(:,:,[100,1100]) = 1;
% 
%     %figure;imagesc(squeeze(artN));
% 
%     Mr(:,nN+1,:) = artN; 

    respVal = mean(Mr(:,:,responseWindow),[1 3]);

    save(sprintf('%s-FFF-respVal',NP.recordingName),'respVal');

    if plotRasters

        bin2= 10;

        M =  BuildBurstMatrix(goodU,round(p.t/bin2),round((stimOn-preBase)/bin2),round((stimDur+preBase*2)/bin2));

        for u =1:size(goodU,2)

            fig = figure;imagesc(squeeze(M(:,u,:)));colormap(flipud(gray(64)));
            xline(preBase/bin2,'k', LineWidth=1.5)
            %Plot stim end:
            xline(stimDur/bin2+preBase/bin2,'k',LineWidth=1.5)
            caxis([0 1])

            cd(NP.recordingDir)
            if ~exist(path+"\Figs",'dir')
                mkdir Figs
            end
            cd(NP.recordingDir + "\Figs")

            print(fig, sprintf('%s-FFF-U%d.png',NP.recordingName,u),'-dpng');

            close 

        end


    end


    if Shuffling_baseline

        if ~isfile(sprintf('FFF-pvalsBaselineBoot-%d-%s.mat',N_bootstrap,NP.recordingName))||repeatShuff==1

            baseline = BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn-preBase)/bin),round((preBase)/bin));
            baseline = single(baseline);
            [nT,nN,nB] = size(baseline);

            % Bootstrapping settings
            N_bootstrap = 1000; % Number of bootstrap iterations
            boot_means = zeros(N_bootstrap, nN,'single');
            resampled_indicesTr = single(randi(nT, [nT, N_bootstrap]));% To store bootstrapped means
            resampled_indicesTi = single(randi(nB, [nB, N_bootstrap]));

            kernel = ones(nT, duration) / (nT * duration); % Normalize for mean
            % Start a parallel pool (if not already started)
            if isempty(gcp('nocreate'))
                parpool; % Start a pool with the default number of workers
            end

            tic
            parfor i = 1:N_bootstrap
                % Resample trials with replacement
                resampled_trials = baseline(resampled_indicesTr(:, i), :,resampled_indicesTi(:, i));
                for ui = 1:nN
                    % Extract the slice for the current unit (t x b matrix)
                    slice = resampled_trials(:, ui, :);
                    slice = squeeze(slice); % Result is t x b

                    % Compute the mean using 2D convolution
                    means = conv2(slice, kernel, 'valid'); % 'valid' ensures the window fits within bounds

                    % Find the maximum mean in this slice
                    boot_means(i, ui) = max(means(:));
                end
            end
            toc


            %%% Calculate p-value & Filter out neurons in which max response window is empty for more than
            %%% 60% of trials

            pvalsResponse = zeros(1,nN);
            ZScoreU = zeros(1,nN);

            %boot_means(:,nN+1) =zeros(1,N_bootstrap);

            for u = 1:nN

                maxWindow = squeeze(Mr(:,u,responseWindow));

                emptyRows = sum(all(maxWindow == 0, 2));

                pvalsResponse(u) = mean(boot_means(:,u)>respVal(u));
                ZScoreU(u) = (respVal(u)-mean(boot_means(:,u)))/(std(boot_means(:,u))+1/(N_bootstrap*nT));

                if emptyRows/nT > 0.6
                    pvalsResponse(u) = 1;
                end

            end
            save(sprintf('FFF-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName),'pvalsResponse')
            save(sprintf('FFF-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName),'ZScoreU')
            save(sprintf('FFF-Base-Boot-%d-%s',N_bootstrap,NP.recordingName),'boot_means')
        end

    end
    %%%%Mean across trials and channels

     if processLFP
         if data.Reference{ex} == "External" && sum(containsFF) >0

             if data.Probe{ex} == "NP1.0"

                 lfpD = NP.getDataLFP([1,100,200,300],stimOn-interStimStats, stimDur+interStimStats*2);

             else

                 lfpD = NP.getData([1,100,200,300],stimOn-interStimStats, stimDur+interStimStats*2);

             end

             [ch,trials,times] = size(lfpD);

             figure;imagesc(squeeze(lfpD(1,:,:)))
             xline((interStimStats/1000)*NP.samplingFrequencyLF)

             for c = 1:ch

                 lfpD(ch,:,:) =  lfpD(ch,:,:)-mean(squeeze(lfpD(ch,:,:)),'all');

             end

         end

     end
     if plotin3D ==1

         if ~exist('figLFP')
             plotBrain
             figLFP = figure(1);
         end
         rLFP = NP.getDataLFP(1:384,stimOn,200);
         rLFP = mean(rLFP,2);
         rLFP = flipud(squeeze(rLFP));
         B = reshape(rLFP, 2, []);
         % Take the mean of every 2 rows
         C = mean(B, 1);
         % Reshape back to a 5x10 matrix
         rLFP = reshape(C, size(rLFP,1)/2, size(rLFP,2));


         cmap = colormap;
         normalizedData = (mean(rLFP,2) - min(mean(rLFP,2))) / (max(mean(rLFP,2)) - min(mean(rLFP,2)));
         colors = interp1(linspace(0, 1, size(cmap, 1)), cmap, normalizedData);

         rLFPred = reshape(mean(reshape(rLFP, 4, []),1),size(rLFP,1)/4,size(rLFP,2));
% 
%          x= linspace(point1(i,1),point2(i,1),size(rLFPred,1));
%          y= linspace(point1(i,2),point2(i,2),size(rLFPred,1));
%          z= linspace(point1(i,3),point2(i,3),size(rLFPred,1));
%          c = mean(rLFPred,2);
         
         %figure;imagesc(mean(rLFPred,2));colorbar

         %Plot units.

         p = NP.convertPhySorting2tIc(NP.recordingDir);

         %Select good units
         label = string(p.label');

         goodU = p.ic(:,label == 'good');
         preBase = 300;

         bin =20;

         [M] = BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn'-preBase)/bin),round((stimDur+2*preBase)/bin));
         Mm = squeeze(mean(M));
         Mm = flipud(Mm);
         fig = figure;imagesc(Mm);colormap(flipud(gray(64)));ylabel('Neurons');xlabel('Time (ms)');
         xline(preBase/bin, '-g', 'Start','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
         xline((preBase/bin+stimDur/bin), '-b', 'End', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
         vectorTicks = 0:100:round(stimDur+preBase*2);
         xticks(vectorTicks/bin);
         xticklabels(vectorTicks-preBase);
         set(gcf,'Color','white');cb = colorbar;
         cb.Label.String = 'Spikes/20 ms';
         
         [Mr]=BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn'+70)/bin),round(50/bin));
         [Mro] = BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn'+stimDur)/bin),round(500/bin));

         res = mean(squeeze(mean(Mr)),2);
         reso = squeeze(mean(Mro));
         figure;imagesc(reso);

         [Mb] =BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn'-preBase)/bin),round(preBase/bin));

         Mb = mean(Mb,3);
         MeanBase = mean(Mb);



     end



    if plotexamplesFF ==1

        rd = NP.getDataLFP(300,stimOn-500,2001);
        
        figure;ylim([min(rd,[],'all') (length(stimOn)+1)*max(rd,[],'all')])
        for s = 1:length(stimOn)
            hold on;plot(s*abs(max(rd,[],'all'))+squeeze(rd(:,s,:)))
        end
        xline([500*NP.samplingFrequencyLF/1000 1500*NP.samplingFrequencyLF/1000])
        vectorTicks = [0 500 1000 1500 2000];
        xticks(vectorTicks*(NP.samplingFrequencyLF/1000))
        xticklabels(vectorTicks-500)
        xlabel('Time (ms)')
        ylabel('Trials (microVolts)')
        prettify_plot

        fig=figure;imagesc(squeeze(rd));c =colorbar;xline(500*NP.samplingFrequencyLF/1000);
        xline([500*NP.samplingFrequencyLF/1000 1500*NP.samplingFrequencyLF/1000])
        vectorTicks = [0 500 1000 1500 2000];
        xticks(vectorTicks*(NP.samplingFrequencyLF/1000))
        xticklabels(vectorTicks-500)
        xlabel('Time (ms)')
        ylabel('Trials')
        c.Ticks = ([-1600:400:800]);
        c.Label.String = 'microVolts';
        set(fig,'Color','white');
        
        
        rd = NP.getDataLFP(1:384,stimOn-500,2001);
        rdm = mean(rd,2);
        verticalDepth = round(unique(sin(deg2rad(data.Angle(ex)))*(data.Depth(ex) - NP.chLayoutPositions(2,:))));
        %take the mean every 2 channels:
        B = reshape(rdm, 2, []);
        % Take the mean of every 2 rows
        C = mean(B, 1);
        % Reshape back to a 5x10 matrix
        rdm = reshape(C, size(rdm,1)/2, size(rdm,3));

        fig=figure;imagesc(flipud(squeeze(rdm)));c =colorbar;xline(500*NP.samplingFrequencyLF/1000);
        xline([500*NP.samplingFrequencyLF/1000 1500*NP.samplingFrequencyLF/1000])
        vectorTicks = [0 500 1000 1500 2000];
        xticks(vectorTicks*(NP.samplingFrequencyLF/1000))
        xticklabels(vectorTicks-500)
        yticklabels(verticalDepth(yticks))
        xlabel('Time (ms)')
        ylabel('Depth (microMeters)')
        %c.Ticks = ([-1600:400:800]);
        c.Label.String = 'microVolts';
        set(fig,'Color','white');
        %set(gca, 'YDir', 'reverse');
        yticks


    end


    %timeSeriesViewer(NP)

    %%% LFP

    %%% Spikes


    %%% shuffling and signifficance

    %%% Depths

    %%% 3D plot

end