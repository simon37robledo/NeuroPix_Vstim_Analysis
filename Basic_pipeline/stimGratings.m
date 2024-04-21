%% PROCESSING STATIC AND DRIFTING GRATINGS STIMULUS


%% Phy commmand: phy template-gui params.py
%Order of batch= FFF, movBall, rectGrid, RectNoise, Gratings.
NumBatch = 1; 

%%

timeSeriesViewer(NP)

redStim = stimOn(1:5:end);

%% Load neuropixels NPAP class and digital TTLs
in =1;

ttl_index = 5; %5 stimulus batches

path = '\\132.66.45.127\data\Large_scale_mapping_NP\Immobilized_exp\PV27\PV27_Experiment_25_6_23\Insertion2\catgt_PV27_Experiment_25_6_23_2_g0';


NP = NPAPRecording(path);


%goodU = p.ic;

[Ttrigger,chNumberT]=NP.getTrigger();

start_end = [Ttrigger{1}];
onset = [Ttrigger{2}];

stimOn = [];
stimOff = [];
stimInter = [];
ttlNum = zeros(1,length(ttl_index));
j =1;

 for i=ttl_index
        stim = onset(onset > start_end(i) & onset < start_end(i+1));

        stimOn = [stimOn stim(1:2:end)]; %general
        stimOff = [stimOff stim(2:2:end)]; %general
        stimInter = [stimInter mean(stimOff-stimOn)];

        ttlNum(j) = length(stim(1:2:end)); %sanity check to see how many stimulus presentations there are per round
        
        j = j+1;
 end

 stimDur = stimOff(1) - stimOn(1);
%% Matlab visual stimulus statistics    
stimDir = '\\132.66.45.127\data\Large_scale_mapping_NP\Immobilized_exp\PV27\PV27_Experiment_25_6_23\Insertion3'; 

file = dir (stimDir);
filenames = {file.name};
gratFiles = filenames(contains(filenames,"StaticDriftingGrating"));
angles = [];
tf = [];
sp = [];

if size(gratFiles) ~= [0 0]
    for i = gratFiles
        grat = load(stimDir+"\"+string(i));
        static_time = cell2mat(grat.VSMetaData.allPropVal(11))*1000; %Static time
        angles = [angles cell2mat(grat.VSMetaData.allPropVal(26))]; %Angles

        tf = [tf cell2mat(grat.VSMetaData.allPropVal(27))]; %time Freq
        sp = [sp cell2mat(grat.VSMetaData.allPropVal(28))]; %spatial Freq
        interStimStats = cell2mat(grat.VSMetaData.allPropVal(40))*1000;

    end
    disp('Visual stats extracted!')
else
    disp('Directory does not exist!');
end

angleNames = unique(angles);
tfNames = unique(tf);
spNames = unique(sp);


%% Load spike sorting and Phy results

% Convert into format to generate raster plots:
p = NP.convertPhySorting2tIc(NP.recordingDir);

%Select good units
label = string(p.label');

goodU = p.ic(:,label == 'good');
amp = p.neuronAmp(label=='good');


disp("Phy extracted!")
%%
cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

%%
%BUILD RASTER PLOT
bin=30;
win=stimDur+interStimStats; 
start = interStimStats/2;
rands = 1000; %Randomization permutations

%image saving path:
pathF = '\\132.66.45.127\data\Large_scale_mapping_NP\Immobilized_exp\PV27\PV27_Experiment_25_6_23\Insertion3\Figures\Gratings';
cd(pathF)

tDiode = "False"; %use diode times?

if tDiode == "True"
   stim = tDiodeOr;
else
    stim = stimOn;
end

    %% 1. BUILD RASTER PLOT: Angles per unit:
    [s, I] = sort(angles);
    angletimesSorted = stimOn(I);

    %Phy name of units:
    c = tsvread(string(NP.recordingDir) + "\cluster_group.tsv");

    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

    for u = 1:length(GoodU_or)


        %Amp = amp(u);
        %nSpikes = goodU(4,u)-goodU(3,u)+1;

        ID = cluster_info.cluster_id(cluster_info.n_spikes==nSpikes & cluster_info.amp == Amp);

        unitCh(u,:) =cluster_info.ch(cluster_info.cluster_id == GoodU_or(u));

        direcNames = unique(directions);

        offsetNames = unique(offsets);

      

        %j = NP.getData(130,500,1000);

        title_stim = sprintf('Insertion-%d-Directions-Unit(%d)-%d-channel-#%d', in, GoodU_or(u),u, cluster_info.ch(cluster_info.cluster_id == GoodU_or(u)));

       
        % Select stimulus
        stim = angletimesSorted - start;
        
        %title_stim = sprintf('Insertion-%d-Angles-Unit-%d-channel-#%s', in, ID, string(goodU(1,u)));

        binR = 30;

        [M]=BuildBurstMatrix(goodU,round(p.t/binR),round(stim/binR),round(win/binR)); %For rasters

        [nTrials,nNeurons,nTimes]=size(M);


       
        %2.STATIC HISTOGRAMS
        trialsA = nTrials/length(angleNames);

        [M2] = BuildBurstMatrix(goodU,round(p.t/bin),round(angletimesSorted/bin),round(static_time/bin)); %For histograms static


        for j = 1:length(angleNames) 
            
            a(j) = sum(squeeze(M2(trialsA*(j-1)+1:trialsA*j,u,:)),"all"); %Plot cricle histogram. 


        end



        %2.1. STATIC: CALCULATE ORIENTATION INDEX
        %Get number of spikes for orthogonal angles to prefered angle

        OIs = sqrt(sum(a.*sin(2*deg2rad(angleNames)))^2 + sum(a.*cos(2*deg2rad(angleNames)))^2)/sum(a); 

        %combine last two angles: (specific for 360 degrees being last index)
        a(1) = (a(1) + a(end))/2;
        a = a(1:end-1);

        
        %2.2. SHUFFLED STATIC
        OIshS = zeros(1,rands);
        aSHs = zeros(1,length(angleNames));


        for s = 1:rands
            
            M_shufS=M2(randperm(size(M2,1)),:,:);
    
            for j = 1:length(angleNames)
    
                aSHs(j) = sum(squeeze(M_shufS(trialsA*(j-1)+1:trialsA*j,u,:)),"all"); %Plot cricle histogram.
    
            end
            
            %2.3. STATIC SHUFFLED: CALCULATE ORIENTATION INDEX
            %Get number of spikes for orthogonal angles to prefered angle
            
   
            %OI formula:
            OIshS(s) = sqrt(sum(aSHs.*sin(2*deg2rad(angleNames)))^2 + sum(aSHs.*cos(2*deg2rad(angleNames)))^2)/sum(aSHs);  % OI calculated as in Dragoi and Swindale. 
            
            %OIshS(s) = (mx + min(aSHs(s,:)) - (OP + OM))/(mx + min(aSHs(s,:)));
    
        end
        
        aSHs = mean(aSHs);
        %combine last two angles: (specific for 360 degrees being last index)
        aSHs(1) = (aSHs(1) + aSHs(end))/2;
        aSHs = aSHs(1:end-1);

        %2.4. STATIC SHUFFLED: PLOT SHUFFLED OI STATISTIC VS LINE WITH REAL
        %OI

        %[h,ps] = ttest(OIshS,OIs, 'Tail','right');

        psR = (sum(OIshS >= OIs) +1)/(rands+1);
        psL = (sum(OIshS <= OIs) +1)/(rands+1);

        if psL < psR
            ps = psL;
        else
            ps = psR;
        end

        figure('visible','off'); 
        if sum(aSHs, 'all') ~= 0
            sh = histfit(OIshS);
            sh(1).FaceColor = "w";
            sh(1).LineStyle= ":";
            title(sprintf('%s -%s trials per angle:',title_stim, string(nTrials/length(angleNames))))
            subtitle('Orientation Index (static).');
            xline(OIs, 'k', LineWidth=1.5);
            text(OIs+0.005,max(sh(1).YData),{'Real OI', sprintf('p = %s', num2str(ps,4))});
            print(sprintf('%s-OI-hist-static.png',title_stim),'-dpng');
            clearvars sh
            close(figure)
            OIstatic(u) = OIs;
        else
            print(sprintf('unit %d does not fire while gratings are static', c(u+1,1)))
        end


        %3. MOVING HISTOGRAMS
        [M3] = BuildBurstMatrix(goodU,round(p.t/bin),round((angletimesSorted + static_time)/bin),round((stimDur-static_time)/bin)); %For histograms moving
        
        for j = 1:length(angleNames) 
            
            a2(j) = sum(squeeze(M3(trialsA*(j-1)+1:trialsA*j,u,:)),"all"); %Plot cricle histogram. 


        end

        %3.1. MOVING: CALCULATE ORIENTATION INDEX
        %Get number of spikes for orthogonal angles to prefered angle


        OIm = sqrt(sum(a2.*sin(2*deg2rad(angleNames)))^2 + sum(a2.*cos(2*deg2rad(angleNames)))^2)/sum(a2); % OI calculated as in Dragoi and Swindale. 

        %combine last two angles: (specific for 360 degrees being last index)

        a2(1) = (a2(1) + a2(end))/2; 
        a2 = a2(1:end-1);

       
        %3.2. SHUFFLED MOVING
        OIshM = zeros(1,rands);
        aSHm = zeros(1,length(angleNames));
        for s = 1:rands
            
            M_shufM=M3(randperm(size(M3,1)),:,:);
    
            for j = 1:length(angleNames)
    
                aSHm(j) = sum(squeeze(M_shufM(trialsA*(j-1)+1:trialsA*j,u,:)),"all"); 
    
            end
            
            %3.3. MOVING SHUFFLED: CALCULATE ORIENTATION INDEX
            %Get number of spikes for orthogonal angles to prefered angle
            

            OIshM(s) = sqrt(sum(aSHm.*sin(2*deg2rad(angleNames)))^2 + sum(aSHm.*cos(2*deg2rad(angleNames)))^2)/sum(aSHm);

            
    
        end
        
        %3.4. MOVING SHUFFLED: PLOT SHUFFLED OI STATISTIC VS LINE WITH REAL

        pmR = (sum(OIshM >= OIm) +1)/(rands+1);
        pmL = (sum(OIshM <= OIm) +1)/(rands+1);

        if pmL < pmR
            pm = pmL;
        else
            pm = pmR;
        end
        
        if sum(aSHm, 'all') ~= 0

            figure('visible','off'); 
            sh = histfit(OIshM);
            sh(1).FaceColor = "w";
            sh(1).LineStyle= ":";
            title(sprintf('%s -%s trials per angle:',title_stim, string(nTrials/length(angleNames))))
            subtitle('Orientation Index (moving).');
            xline(OIm, 'k', LineWidth=1.5);
            text(OIm+0.005,max(sh(1).YData),{'Real OI', sprintf('p = %s', num2str(pm,4))});
            print(sprintf('%s-OI-hist-moving.png',title_stim),'-dpng');
            clearvars sh
            close(figure)
            
            pvals{u} = [ps pm];
        else
            print(sprintf('unit %d does not fire while gratings are moving', ID))
        end

        %4. PLOT RESULTS

        %Plot average accros trials where Units are the y axis. 
        figure('visible','off'); 
        h =imagesc((1:nTimes)*binR,1:2:2*nTrials,squeeze(M(:,u,:)));colormap(flipud(gray(64)));ylabel('Trials');xlabel('Time (ms)');
        title(sprintf('%s -%s trials per angle.',title_stim, string(nTrials/length(angleNames))));
        %cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');clim([0 1]);
        xline(start, '-g', 'Static','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
        xline(start+static_time, '-b', 'Moving', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
        %xline(2000, '-r', 'Inter-trial', 'LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom');
        yline(1:(2*nTrials)/(length(angleNames)):2*nTrials, '--k', string(angleNames), 'LabelVerticalAlignment','bottom', LineWidth=1.5);
        x = 1:2*nTrials/length(angleNames):2*nTrials; %to add first ticks (even);
        y = (nTrials/length(angleNames)):2*nTrials/length(angleNames):2*nTrials; %add second ticks (uneven);
        yticklabels(repmat([1  (nTrials/length(angleNames))/2],1,length(angleNames))); yticks(sort([x y]));
        print(sprintf('%s.png',title_stim),'-dpng');
        clearvars h
        close(figure)

  %Angle histogram per angle (complete stimulus)
%         figure('visible','off'); 
%         for j = 1:length(angleNames) 
% 
%             subplot(length(angleNames),1,j);
%             
%             A(j,:) = sum(squeeze(M(trialsA*(j-1)+1:trialsA*j,u,:))); 
% 
%             
%             %A(j,:) = A(j,:) + (2*nTrials)/(length(angleNames))*j;
% 
%             B = A(j,:)';
% 
%             B = reshape(B,10,round(length(A)/10)); 
% 
%             Af = normalize(sum(B)/mean(sum(B)), 'range');
% 
%             edges = 1:round(win/length(Af)):win;
% 
%             if anynan( Af ) == 0
%                 histogram('BinEdges',edges,'BinCounts',Af, 'FaceAlpha',0.1);
%                 xline(start, 'g', 'S','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='top', LineWidth=2);
%                 xline(start+static_time, '-b', 'M', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='top', LineWidth=2);
%                 ylabel(string(angleNames(j)), 'FontSize',10);
%             end
% 
%         end
% 
%         print(sprintf('%s-hist.png',title_stim),'-dpng');
%         clearvars o
%         close(figure)


%         figure('visible','off'); 
%         j = polarhistogram('BinEdges', deg2rad(angleNames), 'BinCounts', a);
%         title(sprintf('%s -%s trials per angle (static).',title_stim, string(nTrials/length(angleNames))));
%         print(j, sprintf('%s_histS',title_stim),'dpng');
%         clearvars j
%         close(figure)
%         
%         figure('visible','off'); 
%         f = polarhistogram('BinEdges', deg2rad(angleNames), 'BinCounts', a2);
%         title(sprintf('%s -%s trials per angle (moving).',title_stim, string(nTrials/length(angleNames))));
%         print(f, sprintf('%s_histM',title_stim),'dpng');
%         clearvars f
%         close(figure)

    end

    %%
    %2. Time freq per unit:
    [s, I] = sort(tf);
    tfSorted = stimOn(I);

    for u = 1:length(goodU(1,:))
        % Select stimulus
        stim = tfSorted - start;
        title_stim = sprintf('Insertion-%d-Time-freq-Unit-%d-channel-#%s', in, u, string(goodU(1,u)));

        [M]=BuildBurstMatrix(goodU,round(p.t/bin),round(stim/bin),round(win/bin));

        [nTrials,nNeurons,nTimes]=size(M);

        %Plot trials in y axis
        h =imagesc((1:nTimes)*bin,1:2:2*nTrials,squeeze(M(:,u,:)));colormap(flipud(gray(64)));ylabel('Trials');xlabel('Time (ms)');
        title(sprintf('%s -%s trials per time freq.',title_stim, string(nTrials/length(tfNames))));
        %cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');clim([0 1]);
        xline(start, '-g', 'Static','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
        xline(start+static_time, '-b', 'Moving', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
        %xline(2000, '-r', 'Inter-trial', 'LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom');
        yline(1:(2*nTrials)/(length(tfNames)):2*nTrials, '--k', string(tfNames), 'LabelVerticalAlignment','bottom', LineWidth=1.5);
        x = 1:2*nTrials/length(tfNames):2*nTrials; %to add first ticks (even);
        y = nTrials/length(tfNames):2*nTrials/length(tfNames):2*nTrials; %add second ticks (uneven);
        yticklabels(repmat([1  (nTrials/length(tfNames))/2],1,length(tfNames))); yticks(sort([x y]));
        saveas(h,sprintf('%s.png',title_stim));
        clearvars h

    end

    %3. Spatial freq per unit:
    [s, I] = sort(sp);
    spSorted = stimOn(I);

    for u = 1:length(goodU(1,:))
        % Select stimulus
        stim = spSorted - start;
        title_stim = sprintf('Insertion-%d-Spatial-freq-Unit-%d-channel-#%s', in, u, string(goodU(1,u)));

        [M]=BuildBurstMatrix(goodU,round(p.t/bin),round(stim/bin),round(win/bin));

        [nTrials,nNeurons,nTimes]=size(M);

        %Plot trials in y axis
        h =imagesc((1:nTimes)*bin,1:2:2*nTrials,squeeze(M(:,u,:)));colormap(flipud(gray(64)));ylabel('Trials');xlabel('Time (ms)');
        title(sprintf('%s -%s trials per spatial freq.',title_stim, string(nTrials/length(spNames))));
        %cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');clim([0 1]);
        xline(start, '-g', 'Static','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
        xline(start+static_time, '-b', 'Moving', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
        %xline(2000, '-r', 'Inter-trial', 'LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom');
        yline(1:(2*nTrials)/(length(spNames)):2*nTrials, '--k', string(spNames), 'LabelVerticalAlignment','bottom', LineWidth=1.5);
        x = 1:2*nTrials/length(spNames):2*nTrials; %to add first ticks (even);
        y = nTrials/length(spNames):2*nTrials/length(spNames):2*nTrials; %add second ticks (uneven);
        yticklabels(repmat([1  (nTrials/length(spNames))/2],1,length(spNames))); yticks(sort([x y]));
        saveas(h,sprintf('%s.png',title_stim));
        clearvars h

    end


    %4.Summary per neuron (angle):
    
    for a = 1:length(angleNames)
    
        angletimes = stimOn(angles == angleNames(a));

        % Select stimulus
        stim = angletimes - start;
        title_stim = sprintf('Insertion-%d-%s angle', in, string(angleNames(a)));

        
        % Build spike matrix:
        % [M]=BuildBurstMatrix(indexChchannel,t,startTimes,width,class); %t m:

        
        [M]=BuildBurstMatrix(goodU,round(p.t/bin),round(stim/bin),round(win/bin));

        [nTrials,nNeurons,nTimes]=size(M);
        %imagesc(squeeze(M(1,:,:)));colormap(flipud(gray(2)));ylabel('Neurons');xlabel('Time')

        %Plot average accros trials where Units are the y axis. 
        h = imagesc((1:nTimes)*bin,1:2:2*nNeurons,squeeze(mean(M,1)));colormap(flipud(gray(64)));ylabel('Neuron ID');xlabel('Time [ms]');
        title(sprintf('%s -%s trials per angle.',title_stim, string(nTrials)));
        cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');
        xline(start, '-g', 'Static','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom',LineWidth=1.5);
        xline(start+static_time, '-b', 'Moving', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom',LineWidth=1.5);
        yticklabels(string(goodU(1,:))); yticks(1:2:2*nNeurons);
        saveas(h,sprintf('%s.png',title_stim));
        clearvars h
  

    end


    %4.Summary per neuron (time freq.):
    for a = 1:length(tfNames)
    
        tftimes = stimOn(tf == tfNames(a));

        % Select stimulus
        stim = tftimes - start;
        title_stim = sprintf('Insertion-%d-%s-time freq.', in, string(tfNames(a)));

        
        % Build spike matrix:
        % [M]=BuildBurstMatrix(indexChchannel,t,startTimes,width,class); %t m:
  
        [M]=BuildBurstMatrix(goodU,round(p.t/bin),round(stim/bin),round(win/bin));

        [nTrials,nNeurons,nTimes]=size(M);
        %imagesc(squeeze(M(1,:,:)));colormap(flipud(gray(2)));ylabel('Neurons');xlabel('Time')

        %Plot average accros trials where Units are the y axis. 
        h = imagesc((1:nTimes)*bin,1:2:2*nNeurons,squeeze(mean(M,1)));colormap(flipud(gray(64)));ylabel('Neuron ID');xlabel('Time [ms]');
        title(sprintf('%s -%s trials per time freq.',title_stim, string(nTrials)));
        cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');
        xline(start, '-g', 'Static','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
        xline(start+static_time, '-b', 'Moving', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
        yticklabels(string(goodU(1,:))); yticks(1:2:2*nNeurons);
        saveas(h,sprintf('%s.png',title_stim));
        clearvars h
  

    end

    %4.Summary per neuron (spatial freq.):
    for a = 1:length(spNames)
    
        sptimes = stimOn(sp == spNames(a));

        % Select stimulus
        stim = sptimes - start;
        title_stim = sprintf('Insertion-%d-%s-spatial freq.', in, string(spNames(a)));

        
        % Build spike matrix:
        % [M]=BuildBurstMatrix(indexChchannel,t,startTimes,width,class); %t m:

        [M]=BuildBurstMatrix(goodU,round(p.t/bin),round(stim/bin),round(win/bin));

        [nTrials,nNeurons,nTimes]=size(M);
        %imagesc(squeeze(M(1,:,:)));colormap(flipud(gray(2)));ylabel('Neurons');xlabel('Time')

        %Plot average accros trials where Units are the y axis. 
        h = imagesc((1:nTimes)*bin,1:2:2*nNeurons,squeeze(mean(M,1)));colormap(flipud(gray(64)));ylabel('Neuron ID');xlabel('Time [ms]');
        title(sprintf('%s -%s trials per spatial freq.',title_stim, string(nTrials)));
        cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');
        xline(start, '-g', 'Static','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
        xline(start+static_time, '-b', 'Moving', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
        yticklabels(string(goodU(1,:))); yticks(1:2:2*nNeurons);
        saveas(h,sprintf('%s.png',title_stim));
        clearvars h

    end

    %% tuning dependent on depth?
    
    chans = repelem(string(goodU(1,:)),2);
    
    signChans = chans((cell2mat(pvals) <0.05));


%end


