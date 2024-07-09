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
plotin3D=1;

%%
% Iterate through experiments (insertions and animals) in excel file
for ex =15:18%size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
    %ex = 25;
    path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    cd(path)
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
    containsMB = cellfun(@(x) contains(x,'FFF'),Ordered_stims);
    ttlInd = find(containsMB);

    [stimOn stimOff] = NPdiodeExtract(NP,0,"FFF",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff] = NPdiodeExtract(NP,0,"FFF",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A

    stimDur = mean(stimOff-stimOn);
    %%% Raw traces
    %Channel 219 - insertion1 PV67
   
     %%%%Mean across trials and channels

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

         x= linspace(point1(i,1),point2(i,1),size(rLFPred,1));
         y= linspace(point1(i,2),point2(i,2),size(rLFPred,1));
         z= linspace(point1(i,3),point2(i,3),size(rLFPred,1));
         c = mean(rLFPred,2);
         
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