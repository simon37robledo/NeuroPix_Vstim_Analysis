%% stimGratings

cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

bombcellUnits = 0;
newTIC = 0;
repeatShuff =0;
rasters = 1;
shuffling = 0;

newDiode = 0;

bootstraping = 1; 

Shuffling_baseline = 0;
repeatShuff = 0;
tuningIndex =0;
boxplots = 0;
generalResponseSates = 0;

examplesSDG = [8 9 10 11 12 13 14 29 30 31 32 40 41 42 43];
StaticResponses = cell(1,length(examplesSDG));
MovingResponses = cell(1,length(examplesSDG));


%examplesSDG = [8 9 10 11 12 13 14 29 30 31 32];
%%
%SA5_1,PV103_7,PV27_1

k=0;
for ex = 49%1:size(data,1) 7 6 5 40 41 42 43

   
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

%     k=k+1;
%     if ~isempty(StaticResponses{k})
%         %disp()
%         w= sprintf('The following file is already processed: %s',NP.recordingName);
%         warning(w)
%         continue
%     end

    %Create Figs and matData folders if they don't exist

    if ~exist(path+"\Figs",'dir')
        mkdir Figs
    end
% 
%     if ~exist(path+"\matData",'dir')
%         mkdir matData
%     end

 %2. Extract moving ball statistics
    patternIndex = strfind(string(NP.recordingDir), "\catgt");

    endIndex = patternIndex(1)-1;
    stimDir = string(NP.recordingDir);
    stimDir = extractBetween(stimDir,1,endIndex);

    file = dir (stimDir);
    filenames = {file.name};

    file = dir (stimDir);
    filenames = {file.name};
    stimFiles = filenames(contains(filenames,"StaticDrifting"));

    
%         if isempty(stimFiles)
%             %disp()
%             w= sprintf('No static- drifting gratings ball files where found in %s. Skipping into next experiment.',NP.recordingName);
%             warning(w)
%             continue
%         end

    directions = [];
    tempFR = [];
    spatFR = [];

%     stim = load(stimDir+"\"+string(i));
%         static_time = cell2mat(grat.VSMetaData.allPropVal(11))*1000; %Static time
%         angles = [angles cell2mat(grat.VSMetaData.allPropVal(26))]; %Angles
% 
%         tf = [tf cell2mat(grat.VSMetaData.allPropVal(27))]; %time Freq
%         sp = [sp cell2mat(grat.VSMetaData.allPropVal(28))]; %spatial Freq
%         interStimStats = cell2mat(grat.VSMetaData.allPropVal(40))*1000;
    


    if size(stimFiles) ~= [0 0]

        for i = stimFiles %Extract stim properties
            stim= load(stimDir+"\"+string(i));


            directions = [directions cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'angleSequence'))))];

            tempFR = [tempFR cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'tfSequence'))))];

            spatFR = [spatFR cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'sfSequence'))))];

            interStimStats = cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            static_time = cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'static_time'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    Ordered_stims= strsplit(data.VS_ordered{ex},',');
    containsMB = cellfun(@(x) contains(x,'SDG'),Ordered_stims);
    ttlInd = find(containsMB);

    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,newDiode,1,"SDG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,0,1,"SDG",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A

    static_time = cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'static_time'))))*1000; %Static time


    stimInter= mean(stimOn(2:end)-stimOff(1:end-1)); % When dealing with different speeds, save different stim durs, and create M for each speed
    stimDur = mean(-stimOn+stimOff); %When dealing with different speeds, save different stim durs, and create M for each speed

    A = [stimOn directions' tempFR' spatFR'];
    [C indexS] = sortrows(A,[2 3 4]);
    %4. Sort directions:
    directimesSorted = C(:,1)';
    B = [stimOff directions' tempFR' spatFR'];
    [Coff indexSo] = sortrows(B,[2 3 4]);
    directimesSortedOff = Coff(:,1)';

    uDir = unique(directions);nDir = length(uDir);
    uTempFR = unique(tempFR);nTempFR = length(uTempFR);
    uSpatFR = unique(spatFR);nSpatFR = length(uSpatFR);

    if bombcellUnits
        [qMetric,unitType]=NP.getBombCell(NP.recordingDir,1,1);
    end



    %5. Load data to select good units
    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');
    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");
    GoodU_orDepth = cluster_info.depth(cluster_info.group=="good");

    %If new tic matrix needs to be used, move oldTIc files and run convertPhySorting2tIc
    cd(NP.recordingDir)
    if ex == 51
        newTIC =1;
    else
        newTIC =0;
    end

    if isfile("sorting_tIc.mat") && newTIC

        if ~exist('oldTIC', 'dir')
            mkdir oldTIC
        end

        movefile sorting_tIc.mat oldTIC
        movefile sorting_tIc_All.mat oldTIC
    end

    p = NP.convertPhySorting2tIc(NP.recordingDir);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');
    goodSPKS = p.nSpks(label == 'good');
    goodAmp = p.neuronAmp(label == 'good');
    goodUdepth = NP.chLayoutPositions(2,goodU(1,:));

    %6. Get depths of units correcting for angle and micromanipulator depth:
    verticalDepth = sin(deg2rad(data.Angle(ex)))*(data.Depth(ex) - goodUdepth);

    %7. Load raster matrices
    bin = 1;
    preBase = round(stimInter/2);%round(3*interStimStats/4);

    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-stimInter/2)/bin),round((stimDur+stimInter)/bin)); %response matrix
    [nT,nN,nB] = size(Mr);
    trialDivision = nT/(nDir*nSpatFR*nTempFR);

    %figure;imagesc(squeeze(Mr(:,15,:)));colormap(flipud(gray(64)));xline(preBase/bin);xline((preBase+stimDur)/bin);yline(trialDivision:trialDivision:nT-1);

    [MbTotal] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-interStimStats)/bin),round((interStimStats)/bin));% baseline matrix (whole baseline)
    baseline = squeeze(mean(squeeze(mean(MbTotal)),2));
    MrMean = squeeze(mean(Mr));


    stims = stimOn';

    %%%%%%%%%%%%%% Select baseline

    [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round(preBase/bin)); %Baseline matrix plus

    Mbm = mean(Mb,3); %mean across time bins

    spkRateBM = mean(Mbm); %Calculate baseline spike rate of each neuron.

    %2. Calculate denominator for Z-score calculation
    %spkRateBM = mean(Mb); %total mean.
    epsilon = 0.01;
    %denom = mad(MbC,0)+epsilon; %Calculate baseline variation of each neuron.

    denom = mad(Mb,0)+eps;
    direcTrials = trialDivision*nSpatFR*nTempFR;


    trialsPerAngle = nT/nDir;
    
    %%%Response static
    kernel_size = min(500,static_time); %ms

    Ms =  BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted)/bin),round((static_time)/bin));

    statR = zeros(nDir,nN);
    j=1;
    for i =1:direcTrials:nT
        mean_statR = squeeze(mean(Ms(i:direcTrials+i-1,:,100:min(static_time,kernel_size+100)),3));
        statR(j,:) = mean(mean_statR); 
        j=j+1;

    end

    [respValS statPos] = max(statR,[],1);

    %%%Response moving
        %Moving mean approach (works)

    Mm = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted+static_time)/bin),round((stimDur-static_time)/bin));

   
    j = 1;
    movR = zeros(nDir,nN);
    movRpos = zeros(nDir,nN);
    for i =1:trialsPerAngle:nT %iterate through angles
        %         mean_movR = squeeze(mean(movmean(Mm(i:trialsPerAngle+i-1,:,:), kernel_size, 1)));
        %         [movR(j,:), movRpos(j,:)] = max(mean_movR,[],2);
        mean_movR = squeeze(mean(Mm(i:trialsPerAngle+i-1,:,100:min(static_time,kernel_size+100)),3));
        movR(j,:) = mean(mean_movR); 

        j=j+1;
    end

    [respValM, pos] = max(movR,[],1);
    maxMRPos = movRpos(pos);
    %

    if bootstraping
        %%%Bootstrap the baseline per trial
        %%%category and then compare mean b_baseline vs b_response per
        %%%neuron.

        %Bootstrap parameters

        num_trials = direcTrials;
        n_runs = 1000;
        param = 'mean';
        sigmaG = 50;
        duration = 300;
        alpha = 0.05;
        lower_threshold = alpha / 2;
        upper_threshold = 1 - (alpha / 2);
        emptyTrialThres = 0.6;

        statTableM = table();

        statTableS = table();
        tic

        for u = 1:size(goodU,2)


            j = 1;

            baseline = squeeze(Mb(:,u,:));

            bBaseline = get_bootstrapped_equalsamples(baseline,n_runs,num_trials,param);


            for i = 1:direcTrials:nT

                category = uDir(j);

                responseM = squeeze(Mr(i:direcTrials+i-1,u,:));

                responseS = squeeze(Ms(i:direcTrials+i-1,u,:));

                emptyRespM(j) = (sum(sum(responseM,2) > 0))/direcTrials;
                
                emptyRespS(j) = (sum(sum(responseS,2) > 0))/direcTrials;

                

                if sum(bBaseline,'all') ==0 && (sum(responseM,'all')==0 && sum(responseS,'all')==0) %%No spikes in this unit for this stimulus.

                    bootstats_center_base = 0;
                    bootstats_sem_base = 0;
                    activated_Mov = 0;
                    suppresed_Mov = 0;
                    p_boot_Mov = 0;
                    bootstats_center_resp_Mov = 0;
                    bootstats_sem_resp_Mov = 0;
                    activated_Sta = 0;
                    suppresed_Sta = 0;
                    p_boot_Sta = 0;
                    bootstats_center_resp_Sta = 0;
                    bootstats_sem_resp_Sta = 0;
                    validUnit = 0;
                    statTableM = [statTableM; table(u,category,alpha,bootstats_center_base,bootstats_sem_base,activated_Mov,suppresed_Mov,p_boot_Mov,...
                        bootstats_center_resp_Mov,bootstats_sem_resp_Mov)];
                    statTableS = [statTableS; table(activated_Sta,suppresed_Sta,p_boot_Sta,...
                        bootstats_center_resp_Sta,bootstats_sem_resp_Sta, validUnit)];


                else %%% Neuron has spikes during SDG
                    

                    [responseMG]=ConvBurstMatrix(responseM,fspecial('gaussian',[1 duration],sigmaG),'same'); %apply gausian filter to response to smooth spike values

                    % Sham Neuron check
                    %                 responseM = zeros(size(baseline)); %sham neuron
                    %                 baseline = zeros(size(responseM)); %sham neuron
                    %
                    %                 responseM(1:8,1:2:100) = rand(1);
                    %
                    %                 baseline(1:8,1:2:50) = rand(1);

                    [p_boot_Mov, bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(bBaseline,1,responseMG,0,n_runs,num_trials,param);

                    %Calculate probability of bootstats2 >= bootstats1:

                    activated_Mov = p_boot_Mov>upper_threshold;
                    suppresed_Mov = p_boot_Mov<lower_threshold;

                    bootstats_center_base = bootstats_center(1);
                    bootstats_center_resp_Mov = bootstats_center(2);

                    bootstats_sem_base = bootstats_sem(1);
                    bootstats_sem_resp_Mov = bootstats_sem(2);

                    if emptyRespM(j) < emptyTrialThres

                        activated_Mov = 0;

                    end
                

                    statTableM = [statTableM; table(u,category,alpha,bootstats_center_base,bootstats_sem_base,activated_Mov,suppresed_Mov,p_boot_Mov,...
                        bootstats_center_resp_Mov,bootstats_sem_resp_Mov)];


                    %%%%Static


                    [responseSG]=ConvBurstMatrix(responseS,fspecial('gaussian',[1 duration],sigmaG),'same'); %apply gausian filter to response to smooth spike values

                    % Sham Neuron check
                    %                 responseM = zeros(size(baseline)); %sham neuron
                    %                 baseline = zeros(size(responseM)); %sham neuron
                    %
                    %                 responseM(1:8,1:2:100) = rand(1);
                    %
                    %                 baseline(1:8,1:2:50) = rand(1);

                    [p_boot_Sta, bootstats, bootstats_center, bootstats_sem] = get_bootstrap_results_equalsamples(bBaseline,1,responseSG,0,n_runs,num_trials,param);

                    %Calculate probability of bootstats2 >= bootstats1:

                    activated_Sta = p_boot_Sta>upper_threshold;
                    suppresed_Sta = p_boot_Sta<lower_threshold;

                    bootstats_center_base = bootstats_center(1);
                    bootstats_center_resp_Sta = bootstats_center(2);

                    bootstats_sem_base = bootstats_sem(1);
                    bootstats_sem_resp_Sta = bootstats_sem(2);

                    category = uDir(j);

                    validUnit =1;

                    if emptyRespS(j) < emptyTrialThres

                        activated_Sta = 0;

                    end

                    statTableS = [statTableS; table(activated_Sta,suppresed_Sta,p_boot_Sta,...
                        bootstats_center_resp_Sta,bootstats_sem_resp_Sta,validUnit)];

                end


                j = j+1;
            end

        end
        toc

            statTable = [statTableM,statTableS];

           % alpha = 0.05;

            lower_threshold = alpha / 2;
            upper_threshold = 1 - (alpha / 2);


% 
%             statTable.activated_Mov(statTable.activated_Mov ==1) = statTable.p_boot_Mov(statTable.activated_Mov ==1)>upper_threshold;
%             statTable.suppresed_Mov(statTable.suppresed_Mov==1) = statTable.p_boot_Mov(statTable.suppresed_Mov==1)<lower_threshold;
% 
%             statTable.activated_Sta(statTable.activated_Sta == 1) = statTable.p_boot_Sta(statTable.activated_Sta ==1)>upper_threshold;
%             statTable.suppresed_Sta(statTable.suppresed_Sta==1) = statTable.p_boot_Sta(statTable.suppresed_Sta==1)<lower_threshold;

%             for i = 1:direcTrials:nT
% 
%                 responseM = squeeze(Mr(i:direcTrials+i-1,:,:));
% 
%                 responseS = squeeze(Ms(i:direcTrials+i-1,:,:));
% 
%                 emptyRespM = (sum(sum(responseM,3) > 0))/direcTrials;
% 
%                 emptyRespS = (sum(sum(responseS,3) > 0))/direcTrials;
% 
% 
%                 if emptyRespM < emptyTrialThres
% 
%                     statTable.activated_Mov(j) = 0;
% 
%                 end
% 
%                 if emptyRespS(j) < emptyTrialThres
% 
%                     statTable.activated_Sta(j) = 0;
% 
%                 end
% 
%             end


            statSign = statTable(statTable.activated_Sta == 1 | statTable.activated_Mov == 1,:);

            save('bootstrap_stats','statTable');
            save(sprintf('bootstrap_stats_signN-%d',alpha),'statSign')
            RU = unique(statSign.u);
    end

    unique(statSign.u)
%
     if Shuffling_baseline

          N_bootstrap = 1000; % Number of bootstrap iterations
        if ~isfile(sprintf('SDGm-pvalsBaselineBoot-%d-%s.mat',N_bootstrap,NP.recordingName))||repeatShuff==1

            baseline = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-(preBase))/bin),round((preBase)/bin));
            baseline = single(baseline);
            [nT,nN,nB] = size(baseline);

            % Bootstrapping settings
            boot_means = zeros(N_bootstrap, nN,'single');
            resampled_indicesTr = single(randi(nT, [nT, N_bootstrap]));% To store bootstrapped means
            resampled_indicesTi = single(randi(nB, [nB, N_bootstrap]));

            duration = min(300, preBase);

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

                    % Find the maximum mean in this slice (of moving
                    % window)
                    boot_means(i, ui) = max(means(:));
                end
            end
            toc


            %%% Calculate p-value & Filter out neurons in which max response window is empty for more than
            %%% 60% of trials

            pvalsResponseS = zeros(1,nN);
            ZScoreUS = zeros(1,nN);

            pvalsResponseM = zeros(1,nN);
            ZScoreUM = zeros(1,nN);

            %boot_means(:,nN+1) =zeros(1,N_bootstrap);

            for u = 1:nN

                maxWindowS = squeeze(Ms(:,u,:));
                emptyRowsS = sum(all(maxWindowS == 0, 2));

                maxWindowM = squeeze(Mm(:,u,:));
                emptyRowsM = sum(all(maxWindowM == 0, 2));

                pvalsResponseS(u) = mean(boot_means(:,u)>respValS(u));
                ZScoreUS(u) = (respValS(u)-mean(boot_means(:,u)))/(std(boot_means(:,u))+1/(N_bootstrap*trialDivision));

                pvalsResponseM(u) = mean(boot_means(:,u)>respValM(u));
                ZScoreUM(u) = (respValM(u)-mean(boot_means(:,u)))/(std(boot_means(:,u))+1/(N_bootstrap*trialDivision));

                if emptyRowsS/nT > 0.6
                    pvalsResponseS(u) = 1;
                end

                if emptyRowsM/nT > 0.6
                    pvalsResponseM(u) = 1;
                end 

            end
            find(pvalsResponseM<0.05);
            save(sprintf('SDGm-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName),'pvalsResponseM')
            save(sprintf('SDGm-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName),'ZScoreUM')
            save(sprintf('SDGs-pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName),'pvalsResponseS')
            save(sprintf('SDGs-ZscoreBoot-%d-%s',N_bootstrap,NP.recordingName),'ZScoreUS')
            save(sprintf('SDG-Base-Boot-%d-%s',N_bootstrap,NP.recordingName),'boot_means')
        end

    end

    % %Build rasters

    for Rasters = 1
        if rasters ==1
            preR=preBase;
            binr = 50;
            [MrRast] =BuildBurstMatrix(goodU,round(p.t/binr),round((directimesSorted-preR)/binr),round((stimDur+preR*2)/binr));
            %response matrix
            [nT,nN,nB] = size(MrRast);
            
            binr2 = 100;
            [MrRast2] =BuildBurstMatrix(goodU,round(p.t/binr2),round((directimesSorted-preR)/binr2),round((stimDur+preR*2)/binr2));
            %response matrix
            [nT2,nN2,nB2] = size(MrRast2);

            for u = 1:size(goodU,2)
                fig = figure;
                tiledlayout(1,1,"TileSpacing","tight"); %figure('Color','w');
                nexttile
                imagesc(squeeze(MrRast(:,u,:)));colormap(flipud(gray(64)));
                title(sprintf('SDG|%s|UnitN-%d|UnitPhy-%d',strrep(NP.recordingName,'_','-'),u,GoodU_or(u)))
                
              
                xline(preR/binr, '-g', LineWidth=1.5);
                xline(preR/binr+static_time/binr, '-b', LineWidth=1.5);
                xline(nB-preR/binr,'-r',LineWidth=1.5)
                
                ylabel('Trials per angle')
                yline(direcTrials:direcTrials:nT, '--k', LineWidth=1);
                yticks([direcTrials:direcTrials:nT])
                yticklabels([repmat(direcTrials, 1,nDir)]);
                xticks([preR/binr:static_time/binr:nB])
                xticklabels([round(preR/100)*100:static_time:nB*binr])
                %yticklabels(repmat([1  (nTrials/length(tfNames))/2],1,length(tfNames))); yticks(sort([x y]));
                %caxis([0 1])
                yyaxis right
                ylim([1,nT])
                yticks([direcTrials:direcTrials:nT])
                ax = gca;
                ax.YAxis(2).Color = 'k';
                yticklabels(sort(uDir,'descend'))
                ylabel('Angles')
                lims =xlim;
                xt = xticks;
% 
%                 nexttile
% 
%                 y = (mean(squeeze(MrRast2(:,u,:)))./binr2)*1000;
%                 x = 1:nB2; 
%                 y_smooth = smooth(x, y, 0.2, 'loess')';
% 
%                 hold on;
%                 fill([x, fliplr(x)], [y_smooth, zeros(size(y_smooth))], 'b', 'FaceAlpha', 0.3, 'EdgeColor', 'none'); % Fill area
%                 plot(x, y_smooth, '-b', 'LineWidth', 2); % Smoothed curve
%                 plot(x, y, 'ko', 'MarkerFaceColor', 'k'); % Data points
%                 hold off;
%                 xline(preR/binr2, '-g', LineWidth=1.5);
%                 xline(preR/binr2+static_time/binr2, '-b', LineWidth=1.5);
%                 xline(nB2-preR/binr2,'-r',LineWidth=1.5)
%                 xlim([lims(1)*(binr/binr2) lims(2)*(binr/binr2)])
%                 xticks([preR/binr2:static_time/binr2:nB2])
%                 xticklabels([round(preR/100)*100:static_time:nB2*binr2])
%                 try
%                     ylim([min(y)-std(y) max(y)+std(y)])
%                 end
%                 ylabel('Spikes/sec')
%                 xlabel('Time (ms)')
                

% 
%                 h = plot(smoothSum,'LineWidth', 2);
%                 
%                 xticks(round(xt/2))
%                 xticklabels([preR:static_time:nB*binr])
%                 xline(preR/binr/),'-g', LineWidth=1.5);
%                 xline(preR/binr/(size(MrRast,3)/h.NumBins)+static_time/binr/(size(MrRast,3)/h.NumBins)...
%                     , '-b', LineWidth=1.5)
%                 xline(nB/(size(MrRast,3)/h.NumBins)-preR/binr/(size(MrRast,3)/h.NumBins),'-r',LineWidth=1.5)
%                 xlabel('Time (ms)')
%                 yticklabels((round(100*(yticks/size(MrRast,1)))/100)*1000/binr)
%                 ylabel('Spike Rate (spikes/sec)')

                set(gcf,'Color','w')    
                cd(NP.recordingDir+"\Figs")
                pause(1);
                print(fig, sprintf('weirdSDG-%s-Unit-%d.png',NP.recordingName,u),'-dpng');
                clear f
                close 
            end
        end
    end

    if generalResponseSates

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% GENERAL RESPONSE STATISTICS V%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % 1. Statistical test (does the neuron respond to the stimulus in
    % general?)

    %%%% Baseline

    trialsPerAngle = nT/nDir;

    baseR = zeros(nDir,nN);
    baseRpos = zeros(nDir,nN);
    j=1;
    %Moving mean approach (works)
    for i =1:trialsPerAngle:nT %iterate through angles
        mean_baseR = squeeze(mean(movmean(Mr(i:trialsPerAngle+i-1,:,1:round((stimInter/2))/bin), min(kernel_size,stimInter/2), 3)));
        [baseR(j,:), baseRpos(j,:)] = max(mean_baseR,[],2);

        j=j+1;
    end

    [MaxBase, maxBase_position] = max(baseR,[],1);
    maxBase_position = baseRpos(maxBase_position);

    %%%% Static

    %%%Given that neurons have different response delays calculate the min
    %%%or max window across the mean trials applying a convolution.
    kernel_size = 300;
    %Moving kernel convolution approach (didn't work)
%     kernel = ones(1,kernel_size) / kernel_size;
%     meanTR=squeeze(mean(Mr(:,:,round((stimInter/2))/bin:round((stimInter/2))/bin+static_time), 1));
%     
%     % Apply the convolution
%     mean_values = zeros(size(meanTR,1),size(meanTR,2)-length(kernel)+1);
%     for i = 1:size(meanTR,1)%iterate across neurons
%         mean_values(i,:) = conv(meanTR(i,:), kernel, 'valid');
%     end
%     
%     % Find the maximum mean value and its position
%     [staticFR, max_position] = max(mean_values,[],2);
%     
%     % Adjust the position to account for the 'valid' convolution's offset
%     max_position = max_position + floor(kernel_size / 2);
% 
%     %%To calculate appropiate baseline, find max position in baseline
%     %%period (same kernel)
%     meanBase=squeeze(mean(Mr(:,:,1:round((stimInter/2))/bin), 1));
%     mean_valuesBase = zeros(size(meanBase,1),size(meanBase,2)-length(kernel)+1);
%     for i = 1:size(meanTR,1)%iterate across neurons
%         mean_valuesBase(i,:) = conv(meanBase(i,:), kernel, 'valid');
%     end
%     % Find the maximum mean value and its position
%     [MaxBase, maxBase_position] = max(mean_valuesBase,[],2);


    staticFR = zeros(nDir,nN);
    staticFRpos = zeros(nDir,nN);
    j=1;
    %Moving mean approach (works)
    for i =1:trialsPerAngle:nT %iterate through angles
        mean_staticR = squeeze(mean(movmean(Mr(i:trialsPerAngle+i-1,:,round((stimInter/2))/bin:round((stimInter/2))/bin+static_time), kernel_size, 3)));
        [staticFR(j,:), staticFRpos(j,:)] = max(mean_staticR,[],2);

        j=j+1;
    end

    [MaxStaticR, maxSR] = max(staticFR,[],1);
    maxSR = staticFRpos(maxSR);
    observed_diffStatic = MaxStaticR - MaxBase;
   

    %%%% Moving

    
%     %Test of moving window:
%     [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 50],3),'same');
%     figure;u=5;
%     imagesc(squeeze(MrC(:,u,:)));colormap(flipud(gray(64)));xline([maxBase_position(u)-kernel_size/2 maxBase_position(u)+kernel_size/2])
    

    %%Similar to moving ball, scan a window thtough the data at different
    %%orientations, find the max of each orientation, and the total max
    %%across orientations

    movingFRwholeWin = mean(Mr(:,:,round(preBase/bin)+round(static_time/bin):round(preBase/bin)+round(static_time/bin)+round(stimDur/bin)-round(static_time/bin)-round(preBase/bin)),[3 1]);
    observed_diffMovingWW = movingFRwholeWin - spkRateBM;

    movR = zeros(nDir,nN);
    movRpos = zeros(nDir,nN);
    j=1;
    %Moving mean approach (works)
    for i =1:trialsPerAngle:nT %iterate through angles
        mean_movR = squeeze(mean(movmean(Mr(i:trialsPerAngle+i-1,:,...
            round(preBase/bin)+round(static_time/bin):round(preBase/bin)+round(static_time/bin)+round(stimDur/bin)-round(static_time/bin)-round(preBase/bin)), kernel_size, 3)));
        [movR(j,:), movRpos(j,:)] = max(mean_movR,[],2);

        j=j+1;
    end

    [MaxMR, maxMR] = max(movR,[],1);
    maxMR = movRpos(maxMR);
    observed_diffMoving = MaxMR - MaxBase;

    %%%PLOT STATIC R, MOVING R, BASELINE SPKRATES 

    if boxplots
        Responsedata{1} = observed_diffStatic'*1000;
        Responsedata{2} = observed_diffMoving'*1000;
        Responsedata{3} = observed_diffMovingWW'*1000;
        Responsedata{4} = MaxBase'*1000;

        group_names = {'Static response', 'Moving response', 'Moving response WW'};

        h = figure;
        daviolinplot({Responsedata{1,1} Responsedata{1,2}  Responsedata{1,3}} ,'boxcolors','k','outliers',0,'whiskers',1,...
            'box',2,'boxwidth',1.2,'scatter',1,'scattersize',10,'jitter',1,'violinwidth',1.5,'jitterspacing',0.1,...
            'xtlabels', {group_names{1} group_names{2}, group_names{3}});
        title(sprintf('Diff-SDG-%s',strrep(NP.recordingName,'_','-')))
        ylabel('Response - baseline spkRate')
        set(h,'Color','w');%yline(-10:2.5:15,'LineWidth',0.1,'Alpha',0.3);%ylim([-10,15]);
        cd(NP.recordingDir+"\Figs")
        print(h, sprintf('SDG-%s.png',NP.recordingName),'-dpng');
        clear h
    end

    %%% Norm by std of baseline

    stdBaseline = std(squeeze(mean(Mr(:,:,1:round((stimInter/2))/bin),3)));

    zerosSTD = find(stdBaseline==0);
    if ~isempty(zerosSTD)
        stdBaseline(zerosSTD) = std(squeeze(mean(Mr(:,zerosSTD,:),3)));
    end

    StaticResponses{k} = MaxStaticR./stdBaseline;
    MovingResponses{k} = MaxMR./stdBaseline;
   
    %%%
    end

    for Shuffing =1
    if shuffling==1
    rands =1000;
    rand_diffStatic = zeros(rands,nN);
    rand_diffMoving = zeros(rands,nN);

    cd(NP.recordingDir + "\matData")
    if ~isfile(sprintf('SDG-randValues-Static-It-%d.mat',rands)) || repeatShuff ==1

        for s = 1:rands

            %%%Static: Calculate the specific time window (kernel size) where
            %%%the response is max.
            %response
            M_shuff=Mr(:,:,randperm(size(Mr,3)));
            SHstatic = squeeze(mean(M_shuff(:,:,round((stimInter/2)/bin):round((stimInter/2))/bin+round(static_time/bin)),1));

            mean_valuesSS = zeros(size(SHstatic,1),size(SHstatic,2)-length(kernel)+1);
            for i = 1:size(SHstatic,1)%iterate across neurons
                mean_valuesSS(i,:) = conv(SHstatic(i,:), kernel, 'valid');
            end
            SHstatic = max(mean_valuesSS,[],2);

            %baseline
            Mb_shuffStatic = squeeze(mean(M_shuff(:,:,1:round((stimInter/2))/bin), 1));
            mean_valuesB = zeros(size(Mb_shuffStatic,1),size(Mb_shuffStatic,2)-length(kernel)+1);
            for i = 1:size(Mb_shuffStatic,1)%iterate across neurons
                mean_valuesB(i,:) = conv(Mb_shuffStatic(i,:), kernel, 'valid');
            end
            Mb_shuffStatic = max(mean_valuesB,[],2);

            %%%Moving
            SHmoving = mean(M_shuff(:,:,round(preBase/bin)+round(static_time/bin):round(preBase/bin)+round(static_time/bin)+round(stimDur/bin)-round(static_time/bin)-round(preBase/bin)),[3 1]);
            Mb_shuff = mean(M_shuff(:,:,1:preBase),[3 1]);

            rand_diffStatic(s,:) = SHstatic' - Mb_shuffStatic';
            rand_diffMoving(s,:) = SHmoving - Mb_shuff;
        end

        cd(NP.recordingDir + "\matData")
        save(sprintf('SDG-randValues-Static-It-%d',rands),"rand_diffStatic");
        save(sprintf('SDG-randValues-Moving-It-%d',rands),"rand_diffMoving");
    else
        rand_diffStatic = load(sprintf('SDG-randValues-Static-It-%d.mat',rands)).rand_diffStatic;
        rand_diffMoving = load(sprintf('SDG-randValues-Moving-It-%d',rands)).rand_diffMoving;
        rands = size(rand_diffStatic,1);

    end




    p_valueStat = mean(abs(rand_diffStatic) >= abs(observed_diffStatic)');
    p_valueMov = mean(abs(rand_diffMoving) >= abs(observed_diffMoving));



    %%Sort out neurons that respond less than x% of the trials

    exclusionZone = (trialDivision)/nT;
    sumM = sum(Mr,3);
    exludedN = [];

    inclusionZone = trialDivision*0.7; %Number of trials in the same trial category that are enough as a good unit.

    for u =1:nN

        if sum(sumM(:,u) > 0)<floor(nT*exclusionZone)

            exludedT = ones(1,nT/trialDivision);

            for i=1:trialDivision:nT

                if sum(sumM(i:i+trialDivision-1,u)) <= inclusionZone

                    exludedT(i) = 0;
                end

            end

            if sum(exludedT) >= 1

                exludedN = [exludedN u];

            end

        end

    end

    %%% plot rand diff values as histogram, then real values as lines.
    %%% Then at each side plot a red line where signifficance is
    %%% crossed.

    %1. normalize across neurons and iterations random events:

    NstR = normalize(normalize([rand_diffStatic;observed_diffStatic'],'center','scale',1),2);

    NmvR = normalize(normalize([rand_diffMoving;observed_diffMoving],'center','scale',1),2);

    %2. normalize observed values across neurons:

    Nst = normalize(observed_diffStatic);
    Nmv = normalize(observed_diffMoving);

    %%ViolinPlot

    Responsedata{1} = NstR(end,:)';
    Responsedata{2} = mean(NstR(1:end-1,:),1)';
    Responsedata{3} = NmvR(end,:)';
    Responsedata{4} = mean(NmvR(1:end-1,:),1)';

    group_names = {'Obs. SR', 'Rand. SR' , 'Obs. MR', 'Rand. MR'};

    figure;
    h = daviolinplot({Responsedata{1,1} Responsedata{1,3}} ,'boxcolors','k','outliers',0,'whiskers',1,...
    'box',2,'boxwidth',0.8,'scatter',1,'scattersize',10,'jitter',1, 'jitterspacing',0.1,...
    'xtlabels', {group_names{1} group_names{3}});


    %%%PVALUE HISTOGRAM

    se_valuesST = std(rand_diffStatic, 0, 2) / sqrt(size(rand_diffStatic, 2));
    se_valuesMV = std(rand_diffMoving, 0, 2) / sqrt(size(rand_diffMoving, 2));

    rh=histogram(mean(NstR(1:end-1,:),2),'FaceColor','w','EdgeColor','k');
    edges = [min(rh.BinEdges)-mean(diff(rh.BinEdges))*2,min(rh.BinEdges)-mean(diff(rh.BinEdges)), rh.BinEdges,...
        max(rh.BinEdges)+mean(diff(rh.BinEdges)), max(rh.BinEdges)+mean(diff(rh.BinEdges))*2];
    %create plotable observed values for static and moving
    %Static:
    plotObs = NstR(end,:);
    plotObs(plotObs> max(rh.BinEdges)) = max(rh.BinEdges)+mean(diff(rh.BinEdges))/2;
    plotObs(plotObs<min(rh.BinEdges)) = min(rh.BinEdges)-mean(diff(rh.BinEdges))/2;
    %Moving:
    plotObsm = NmvR(end,:);
    plotObsm(plotObsm> prctile(mean(NstR(1:end-1,:),2), 99.95)) = max(rh.BinEdges)+mean(diff(rh.BinEdges)*1.5);
    plotObsm(plotObsm<prctile(mean(NstR(1:end-1,:),2), 0.05)) = min(rh.BinEdges)-mean(diff(rh.BinEdges)*1.5);
    close gcf


    figure('Color','w'); %%%%%To do: Why I don't find real responses as within thwe 99 random percentile??
    %plot histogram with random responses
    yyaxis right;histogram(mean(NstR(1:end-1,:),2));hold on;histogram(mean(NmvR(1:end-1,:),2),FaceAlpha=0.2);
    xline(prctile(mean(NstR(1:end-1,:),2), 99.95),'k','Upper 0.05 tail','LabelHorizontalAlignment', 'left',LabelVerticalAlignment='middle');
    xline(prctile(mean(NstR(1:end-1,:),2), 0.05),'k','Lower 0.05 tail','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='middle');
    ylabel('Randomized neurons')
    xlabel('(Response - Baseline)')

    %Plot histogram with observed values in left axis
    yyaxis left;histogram(plotObs,'BinEdges',edges);hold on; histogram(plotObsm,'BinEdges',edges,FaceAlpha=0.2);
    ylabel('Observed neurons')

    legend(string(sprintf('%d Ob. Neur. SR.',nN)),string(sprintf('%d Ob. Neur. MR.',rands)),...
        string(sprintf('%d iter. SR.',rands)),string(sprintf('%d iter. MR.',rands)));

    hold on; histogram(plotObsm(plotObsm< max(rh.BinEdges)),'BinEdges',edges)

    plotObsm(plotObsm>min(prctile(mean(NstR(1:end-1,:),2), 0.05)) & plotObsm<max(prctile(mean(NstR(1:end-1,:),2), 99.95)))

    %create plotable observed values

    %         pltObs = NstR;
    %         plotObs(plotObs>prctile(mean(NstR(1:end-1,:),2), 99.5) ) =
    %
    %         && plotObs>xline(prctile(mean(NstR(1:end-1,:),2), 0.5)
    %         xline(NstR(end,:))
    %Produce hist with data such that if the value is above the
    %percentiles, it has a fixed value.
    modObsST = observed_diffStatic;
    %modObsST(observed_diffStatic>) = observed_diffStatic()


    hold on; xline(log10(observed_diffStatic+abs(min(observed_diffStatic))+1));
    xline(Nst)
    hold on;
    %
    %         % Add error bars for the standard error
    %         errorbar(1:length(mean(se_valuesST,2)), mean(se_valuesST,2), se_valuesST, 'k.', 'LineWidth', 1);

    end
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%% TUNING INDEX STATISTICS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % 3. MOVING SECTION

    for tuning =1

        if tuningIndex ==1
            uDir = unique(directions);
            uDirRad = deg2rad(uDir);

            [M3] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted + static_time)/bin),round((stimDur-static_time)/bin)); %For histograms moving


            trialDivision = nT/(length(uDir));

            aM = zeros(nN,length(uDir));
            for u =1:nN
                for j = 1:length(uDir)
                    aM(u,j) = sum(squeeze(M3(trialDivision*(j-1)+1:trialDivision*j,u,:)),"all"); %Plot cricle histogram.
                end
            end


            OIm = sqrt(sum(aM.*sin(2*uDirRad),2).^2 + sum(aM.*cos(2*uDirRad),2).^2)./sum(aM,2);


            %3.2. SHUFFLED MOVING
            rands =1000;
            OIshM = zeros(nN,rands);
            aSHm = zeros(nN,length(uDir));
            for s = 1:rands

                M_shufM=M3(randperm(size(M3,1)),:,:);

                trialDivision = nT/(length(uDir));

                for u =1:nN
                    for j = 1:length(uDir)
                        aSHm(u,j) = sum(squeeze(M_shufM(trialDivision*(j-1)+1:trialDivision*j,u,:)),"all"); %Plot cricle histogram.
                    end
                end

                OIshM(:,s) = sqrt(sum(aSHm.*sin(2*uDirRad),2).^2 + sum(aSHm.*cos(2*uDirRad),2).^2)./sum(aSHm,2);

            end

            %3.4. MOVING SHUFFLED: PLOT SHUFFLED OI STATISTIC VS LINE WITH REAL

            pmR = (sum(OIshM >= OIm,2) +1)/(rands+1);
            pmL = (sum(OIshM <= OIm) +1)/(rands+1);

            %%Rand
            rand = 1000;

            bin2 =50;

            %%%%Select good responsive neurons:





            %%%%%%

            [Mr2] = BuildBurstMatrix(goodU,round(p.t/bin2),round((directimesSorted-stimInter/2)/bin2),round((stimDur+stimInter)/bin2)); %response matrix
            figure;
            imagesc(squeeze(Mr2(:,1,:)));colormap(flipud(gray(64)));
            %Plot stim start:
            xline(preBase/bin2,'k', LineWidth=1.5)
            %Plot stim end:
            xline(stimDur/bin2+preBase/bin2,'k',LineWidth=1.5)
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

        end
    end



end

%% Summary statists

% for i = 1:7
%     
% 
% end
% 
% responseStat{1} = meanIn';
% for i = 1:7
%    
% 
% end
% responseStat{2} = meanIn';

%group_names = {'staticR' 'movingR'};

animalColor = [repmat([0.4940 0.1840 0.5560],7,1);repmat([0.4660 0.6740 0.1880],8,1);...
    repmat([0.3010 0.7450 0.9330],4,1);repmat([0.6350 0.0780 0.1840],4,1)];

figure
hold on

animalNames = {'PV103','PV27','SA5','PV35'};

j=2;
yt =3;
name =1;
for i = 1:length(examplesSDG)

    j=j+1;

    SR = StaticResponses{i};
    MR = MovingResponses{i};

    indTopSR = find(SR>=prctile(SR,75));
    meanInS = mean(SR(SR>=prctile(SR,75)),'omitnan');
    meanInM = mean(MR(indTopSR),'omitnan');

    plot([1,2],[meanInS,meanInM],'-o','MarkerFaceColor',animalColor(i,:),'MarkerEdgeColor',animalColor(i,:),'Color',animalColor(i,:));

    if animalColor(min(j,length(animalColor)),1) ~= animalColor(j-1,1) || j>length(animalColor)
        yt=yt+4;
        text(2.5, yt, animalNames{name}, 'Color', animalColor(i,:), 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
         name = name+1;
    end

end
hold off
xlim([0.5,2.5])
xlabel('static Response VS Moving response')
ylabel('mean Normalized response across neurons (std)')
title('Static vs Moving, Gratings. Line = insertion, Color = animal')

%%
%%%Moving
j=2;
yt =3;
name =1;

figure
hold on
for i = 1:length(examplesSDG)

    j=j+1;

    SR = StaticResponses{i};
    MR = MovingResponses{i};

    indTopMR = find(MR>=prctile(MR,75));
    meanInM = mean(MR(MR>=prctile(MR,75)),'omitnan');
    meanInS = mean(SR(indTopMR),'omitnan');

    plot([1,2],[meanInS,meanInM],'-o','MarkerFaceColor',animalColor(i,:),'MarkerEdgeColor',animalColor(i,:),'Color',animalColor(i,:));

    if animalColor(min(j,length(animalColor)),1) ~= animalColor(j-1,1) || j>length(animalColor)
        yt=yt+4;
        text(2.5, yt, animalNames{name}, 'Color', animalColor(i,:), 'FontSize', 12, ...
        'FontWeight', 'bold', 'HorizontalAlignment', 'left', 'VerticalAlignment', 'bottom');
         name = name+1;
    end

end
hold off
xlim([0.5,2.5])
xlabel('static Response VS Moving response')
ylabel('mean Normalized response across neurons (std)')
title('Static vs Moving, Gratings. Line = insertion, Color = animal')
% Add the legend (only for the color groups)



%legend(animalNames,"Color",mat2cell(unique(animalColor,'rows'),ones(1,size(unique(animalColor,'rows'),1)),3))




