%% Overall response score

basic_pathPV102 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV102';
expPV102 = 'PV102_experiment_18_7_23';

basic_pathPV103 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV103';
expPV103 = 'PV103_Experiment_12_6_23';

basic_pathPV67 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV67';
expPV67= 'PV67_experiment_5_7_23';

basic_pathPV27 = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV27';
expPVPV27 = 'PV27_Experiment_25_6_23';

basicPathA = {basic_pathPV67,basic_pathPV103,basic_pathPV27};

expA = {expPV67,expPV103,expPVPV27};

StimOrder = {'fullFieldFlash','linearlyMovingBall','rectGrid','StaticDriftingGrating'};

sPV67 = [1,2,3 ;find(strcmp(StimOrder, 'fullFieldFlash')==1),find(strcmp(StimOrder, 'linearlyMovingBall')==1),find(strcmp(StimOrder, 'rectGrid')==1)];
sPV103 = [1,2,3;find(strcmp(StimOrder, 'rectGrid')==1),find(strcmp(StimOrder, 'linearlyMovingBall')==1),find(strcmp(StimOrder, 'StaticDriftingGrating')==1)];
sPV27 = [1,2,3;find(strcmp(StimOrder, 'rectGrid')==1),find(strcmp(StimOrder, 'linearlyMovingBall')==1),find(strcmp(StimOrder, 'StaticDriftingGrating')==1)];

plotexamples = 0;
plotexamplesRG =0;

stimTTL = {sPV67,sPV103,sPV27};


MinusRBt = cell(4,length(expA));

MinVal = repmat(inf,1,length(StimOrder));
MaxVal = repmat(-inf,1,length(StimOrder));

% Vals = cell(3,length(StimOrder));
% Vals(1,:) = StimOrder;
% Vals
bin = 50; %response bin ms

preBase = 500;

%%
preBase2 = 500;

binsPB = preBase2/bin; %Baseline bins.


MrSDG = zeros(18,100+binsPB*2); %4 bin*4 for prebase

MrLMB = zeros(18,47+binsPB*2);

MrFFF = zeros(18,21+binsPB*2);

MrRG = zeros(18,21+binsPB*2);

binsBase = preBase/bin;

%%
dbstop if warning

itN =1;

for a = 1:length(expA)

     MinusRBi = cell(1,length(animalA{a}));
     Zscorei = cell(1,length(animalA{a}));
     Tunei = cell(1,length(animalA{a}));

    for in = 1:length(animalA{a})

       

         if a == 3
                
            path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
            (in+1)+string(filesep)+"catgt_"+string(expA{a})+"_"+(in+1)+"_g0");
        else

            path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
                in+string(filesep)+"catgt_"+string(expA{a})+"_"+in+"_g0");
        end


        NP = NPAPRecording(path);

        cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

        %Good units

        GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

        p = NP.convertPhySorting2tIc(NP.recordingDir);

        %Select good units
        label = string(p.label');

        goodU = p.ic(:,label == 'good');

        if length(GoodU_or) ~= length(goodU)

            disp(NP.recordingDir)
        end


        %Stim times

        [Ttrigger,chNumberT]=NP.getTrigger();

        Sstart = [Ttrigger{1}];
        Send = [Ttrigger{2}];
        onset = [Ttrigger{3}];
        offset = [Ttrigger{4}];

        MinusRBs = cell(2,length(stimTTL{a}(1,:)));

        Zscores = cell(2,length(stimTTL{a}(1,:)));

        tune = cell(2,length(stimTTL{a}(1,:)));

        %%%Example neurons:

        if a == 1 & in == 1

            eNeuron = 115:125;

        elseif a == 2 & in == 1

            eNeuron = 2:12;

        else a == 3 & in ==1

            eNeuron = 98:108;

        end

        for s = 1:length(stimTTL{a}(1,:))

            if a == 3 & length(Send)>3
                Send = Send(2:end); %PV27 weirdness, it has one more ttl off trigger at the start
            end
    

            ttl_index = stimTTL{a}(1,s);

            stimName  = StimOrder{stimTTL{a}(2,s)};

            stimOn = [];
            stimOff = [];
            stimInter = [];
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


            %%%%%%%Baseline

            stims = stimOn;

            %Select biggest response baseline
            
            basePresp = preBase+stimDur;

            [Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round((preBase)/bin));

            Mb = mean(Mb,3);

            %find(GoodU_or == 53)
%             figure
%             imagesc(squeeze(MbC(:,120)));colormap('Gray')
%             colorbar
%             xlabel('1 bin (500 ms)')
%             ylabel('trials')
%             title('Baseline spike rate of example neuron mov ball, insertion 1 experiment 1')

            %get stiminter starting at bin ms after the previous stim.

            [trials, neurons] = size(Mb);

%             MbC = zeros(trials, neurons);

            percentTrials = round(trials*0.2);

%             % Apply medfilt1 along the first dimension (rows)
%             for i = 1:neurons
%                 MbC(:, i) = medfilt1(Mb(:, i), percentTrials);
%             end
           

            [S, Sindex] = sort(Mb);         

           % Mb = S(1:end-percentTrials,:);

            [trials, neurons] = size(Mb);

            %randomize again.
            Mb = Mb(randi([1, size(Mb,1)], 1, size(Mb,1)),:);

            %Take mean across window of % of trials

            MbC = zeros(round(trials/percentTrials), neurons); 

            for i = 1:percentTrials:trials


                meanb = mean(Mb(i:min(i+percentTrials-1, end),:),1);

                MbC(j,:) = meanb;

                j = j+1;

            end


           % Mb = Mb(Sindex);

            spkRateBM = mean(MbC);

%             figure
%             imagesc(squeeze(Mb(:,50)));colormap('Gray')
%             colorbar
%             xlabel('1 bin (500 ms)')
%             ylabel('trials')
%             title('Baseline spike rate of responsive neuron insertion 2')
            
            %%%%Response

            %Select biggest response baseline

            [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims)/bin),round((stimDur)/bin)); %get stimdur plus bin ms after offset.

            [nT,nN,nB] = size(Mr);
  
            
            %%%%% Select time bin of response per stimulus:
% 
           
               
       %%%%%% Static and drifting grating:  
            if length(stimName) == length('StaticDriftingGrating')

                  % Matlab visual stimulus statistics
                if a == 3 %PV27

                    stimDir =convertStringsToChars(string(basicPathA{a})+"\"+string(expA{a})+"\Insertion"+(in+1));
                else

                    stimDir =convertStringsToChars(string(basicPathA{a})+"\"+string(expA{a})+"\Insertion"+in);
                end

                file = dir (stimDir);
                filenames = {file.name};
                recFiles = filenames(contains(filenames,"StaticDriftingGrating"));
                angles = [];
                tf = [];
                sp = [];


                if size(recFiles) ~= [0 0]

                    j =1;
                    for i = recFiles
                        grat = load(stimDir+"\"+string(i));
                        static_time = cell2mat(grat.VSMetaData.allPropVal(11))*1000; %Static time
                        angles = [angles cell2mat(grat.VSMetaData.allPropVal(26))]; %Angles

                        tf = [tf cell2mat(grat.VSMetaData.allPropVal(27))]; %time Freq
                        sp = [sp cell2mat(grat.VSMetaData.allPropVal(28))]; %spatial Freq
                        interStimStats = cell2mat(grat.VSMetaData.allPropVal(40))*1000;

                        j = j+1;
                    end
                    disp('Visual stats extracted!')
                else
                    disp('Directory does not exist!');
                end

                angleNames = unique(angles);

               
                [sA, I] = sort(angles);
                
                stims = stims(I);

                [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round((stimDur+preBase*2)/bin)); %get stimdur plus bin ms after offset.

                [nT,nN,nB] = size(Mr);


                binselectStart = 200/bin; %200 ms

                binselectEnd = 700/bin; %700 ms
 

                spkRateR = mean(Mr(:,:,binsBase+binselectStart:binsBase+binselectEnd),[1 3]);
                

                epsilon = 0.01;

                denom = mad(MbC,0)+epsilon; %mean(Mb,0)+epsilon; %

                mSpk = mean(spkRateR);

                Zscore = (spkRateR - (spkRateBM + spkRateR)/2)./denom;

                max(max_position)%./(spkRateR+spkRateB);

                sDirections = sort(directions);

               % tunning = sDirections(max_position(:,1));

                %MinusRB = (mean(spkRateR) -(mean(spkRateB)))./(mean(spkRateR) + (mean(spkRateB)));

                MinusRB = spkRateR - spkRateBM;

                plotExamplesSDG =0;


                if (a == 1 && in ==1) || (a == 2 && in ==1) || (a == 3 && in ==1)


                    if plotExamplesSDG ==1

                        for u = eNeuron

                            figure%('visible','off');
                            h =imagesc((1:nT)*bin,1:2:2*nT,squeeze(Mr(:,u,:)));colormap(flipud(gray(64)));ylabel('Trials');xlabel('Time (ms)');
                            title(sprintf('In.%d-Gratings-raster-U%d-S%.3f',in,u,Zscore(u)))
                            %cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');clim([0 1]);
                            xline(preBase, '-g', 'Static','LabelHorizontalAlignment', 'left',LabelVerticalAlignment='bottom', LineWidth=1.5);
                            xline(preBase+static_time, '-b', 'Moving', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
                            xline([preBase+4*bin,preBase+14*bin]);
                            %xline(2000, '-r', 'Inter-trial', 'LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom');
                            yline(1:(2*nT)/(length(angleNames)):2*nT, '--k') %string(angleNames), 'LabelVerticalAlignment','bottom', LineWidth=1.5);
                            xline(round(preBase+stimDur),LineWidth=1.5)
                            x = 1:2*nT/length(angleNames):2*nT; %to add first ticks (even);
                            y = (nT/length(angleNames)):2*nT/length(angleNames):2*nT; %add second ticks (uneven);
                            yticklabels(repmat([1  (nT/length(angleNames))/2],1,length(angleNames))); yticks(sort([x y]));
                            prettify_plot
                            print(sprintf('In.%d-Gratings-raster-U%d.png',in,u),'-dpng');
                            clearvars h
                            close
                        end

                    end
                end

                        
       %%%%%% Full field flash:  
            elseif length(stimName) == length('fullFieldFlash')

                tunning = 0;
                
                spkRateR = mean(Mr(:,:,2:12),[1 3]);

                epsilon = 0.01;

                denom = mad(MbC,0)+epsilon; %mean(Mb,0)+epsilon; %

                mSpk = mean(spkRateR);

                Zscore = (spkRateR - (spkRateBM + spkRateR)/2)./denom;
                          

                %tunning = sDirections(max_position(:,1)); 


                %MinusRB = (mean(spkRateR) -(mean(spkRateB)))./(mean(spkRateR) + (mean(spkRateB)));

                MinusRB = spkRateR - spkRateBM;


       %%%%%% Ractangle grid:        
            elseif length(stimName) == length('rectGrid')


                if a ==1 && in == 1
                    2+2

                end
                

                % Matlab visual stimulus statistics
                if a == 3

                    stimDir =convertStringsToChars(string(basicPathA{a})+"\"+string(expA{a})+"\Insertion"+(in+1));
                else

                    stimDir =convertStringsToChars(string(basicPathA{a})+"\"+string(expA{a})+"\Insertion"+in);
                end

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

                [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((stims-preBase)/bin),round((stimDur+preBase*2)/bin)); %get stimdur plus bin ms after offset.

                [nT,nN,nB] = size(Mr);

                [~,indRG] = sort(seqMatrix);

                trialsPerCath = length(seqMatrix)/(length(unique(seqMatrix)));

                %Take mean across cathegories

                MaxResp = zeros(nT/(trialsPerCath), nN);

                ZscorePos = zeros(nT/(trialsPerCath), nN);

                %sort Mr accordyinh to stimulus type (position)
                Mr = Mr(indRG,:,:);

                MrC = ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');

                  %%%% Remove outliers to cut noise neurons

                  binsBase = preBase/bin;

                    j = 1;

                    for i = 1:trialsPerCath:nT

                        OnePosMatrix = MrC(i:min(i+trialsPerCath-1, end),:,:);

 
                        
                        %%% On response
                        meanBinON = mean(OnePosMatrix(:,:,binsBase+2:binsBase+8),3); %onset response

                        meanTrialON = mean(meanBinON,1);

                        stdON = std(meanBinON,1);

                       
                        %%% Off response
                        meanBinOFF = mean(OnePosMatrix(:,:,binsBase+(round((stimDur)/bin)+1):end-2),3); %offset response

                        meanTrialOFF = mean(meanBinOFF,1);

                        stdOFF = std(meanBinOFF,1);

                        %%% Baseline mean
                        meanBinBase = mean(OnePosMatrix(:,:,1:binsBase),3);

                        meanTrialBase = mean(meanBinBase,1);

                        stdBase = std(meanBinBase,1);



                        for u = 1:nN

%                             if i == 1 & u == 214
%                                 2+2
%                             end

                            %%% Calculate outliers for On, Off, and
                            %%% Baseline. 
                            
                            OutsOn = meanBinON(:,u) < meanTrialON(u)+2*stdON(u); %Remove outlier trials
                            NewMeanON = mean(meanBinON(OutsOn,u));
                            NewMeanON(isnan(NewMeanON)) = 0;

                            OutsOFF = meanBinOFF(:,u) < meanTrialOFF(u)+2*stdOFF(u);
                            NewMeanOFF = mean(meanBinOFF(OutsOFF,u));
                            NewMeanOFF(isnan(NewMeanOFF)) = 0;

                            OutsBase = meanBinBase(:,u) < meanTrialBase(u)+2*stdBase(u);
                            NewMeanBase = mean(meanBinBase(OutsBase,u));
                            NewMeanBase(isnan(NewMeanBase)) = 0;

                            [MaxRespi, OffOrOn]= max([NewMeanOFF,NewMeanON]); %Select max between onset and offset response. 
                            
                            MaxResp(j,u) = MaxRespi;

                            epsilon = 0.02;

                            ZscorePos(j,u) = (MaxResp(j,u)-NewMeanBase)/(stdBase(u)+epsilon);

%                             figure
%                             imagesc(squeeze(OnePosMatrix(:,u,:)))
%                             colorbar


                        end

                        j = j+1;

                    end

                    denom = mad(MbC,0)+epsilon;

                    [~, maxSZIndex]= max(ZscorePos);

                    tunning = positionsMatrix(maxSZIndex,:)./max(positionsMatrix(:));


%                     for u =1:nN
%                         if sum(ZscorePos(:,u)>2,1)>1
% 
%                             [pos, uni] = find(ZscorePos(:,u)>2);
% 
%                             tunning(u,:) = mean(positionsMatrix(pos,:),1)./max(positionsMatrix(:));
% 
%                         end
%                     end

                    [MaxZS,maxSZIndex]= max(ZscorePos);

                    Zscore = MaxZS;


                    plotexamplesRG = 0;

                    if a == 3 && in ==2
                        2+2
                    end

                    if (a == 1 && in ==1) || (a == 2 && in ==1) || (a == 3 && in ==2)


                        if plotexamplesRG == 1


                           eNeuron = 1:nN;

                           eNeuron = 195:203;

                            for u = eNeuron

                                

                                %%%%%%Heat map
                                ZscoreScreen = reshape(ZscorePos(:,u),max(positionsMatrix));

                                %ZscoreScreen(ZscoreScreen<3) = 0;
 

                                [sx sy] = find(ZscoreScreen>2);


                                fig = figure;
                                imagesc(ZscoreScreen);colormap(flipud(gray));
                                c = colorbar;
                                caxis([0 5]);
                                c.Title.String = "Z-score";
                                xticks([1:length(ZscoreScreen)])
                                yticks([1:length(ZscoreScreen)])
                                hold on
                                plot(sy,sx,'r.','MarkerSize',20)
                                xline(1.5:1:length(ZscoreScreen)+0.5)
                                yline(1.5:1:length(ZscoreScreen)+0.5)
                                ylabel('Y screen pos.')
                                xlabel('X screen pos.')
                                prettify_plot
                                print(fig, sprintf('%s-In.%d-rect-GRid-dir-U%d.png',expA{a},in,u),'-dpng')
                                close

                                %%%%%Raster

                                t = tiledlayout(sqrt(max(seqMatrix)), sqrt(max(seqMatrix)),'TileSpacing','tight');


                                PositionsTotal = positionsMatrix(sort(seqMatrix),:);
                                [posS,indexX] = sortrows(PositionsTotal,1); %Sort first dimension because tile layout moves through columns

                                MrS = Mr(indexX,:,:);
                                
                                for i = 1:trialsPerCath:nT
                                    %Build raster
                                    M = MrS(i:min(i+trialsPerCath-1, end),u,:);
                                    [nTrials,nNeurons,nTimes]=size(M);
                                    nexttile
                                    imagesc((1:nTimes)*bin,1:nTrials,squeeze(M));colormap(flipud(gray(64)));
                                    xline(preBase, '--', LineWidth=2, Color="#77AC30");
                                    xline(stimDur+preBase, '--', LineWidth=2, Color="#0072BD");
                                    xticks([preBase round(round(stimDur/100))*100+preBase]);
                                    caxis([0 1]);
                                    set(gca,'YTickLabel',[]);
                                    
                                    if i < nT - trialsPerCath*max(positionsMatrix(:))-1
                                        set(gca,'XTickLabel',[]);

                                    end
                                end
                                fig = gcf;
                                set(fig, 'Color', 'w');
                                colorbar
                                % Set the color of the figure and axes to black
                                title(t,sprintf('In.%d-rect-GRid-raster-U%d',in,u))
                                ylabel(t,sprintf('%d trials',nTrials))
                                fig.Position = figPosRG;
                                prettify_plot
                                print(sprintf('%s-In.%d-rect-GRid-raster-U%d.png',expA{a},in,u),'-dpng')
                                close

                            end

                        end
                    end


                    [MaxRespS maxIndex]= sort(MaxResp,'descend'); %% sort response strengt

                    spkRateR = MaxRespS(1,:);

                    MinusRB = spkRateR - spkRateBM;








               % [] = find(HM ==83)

               


                %%%% For responsive neurons plot receptive field

                Si = any(ZscorePos > 3);

                ZscorePos(ZscorePos>5) = 5;

                ZscorePos(ZscorePos<3) = 0;
% 
%                 figure;
%                 imagesc(ZscorePos(:,Si))
%                 colorbar
% 
%                 mean(MaxRespS(:,Si(j)))
%                 xlabel('neuron')
%                 ylabel('Screen position')

                
            
       %%%%%% Moving ball:    
            elseif length(stimName) == length('linearlyMovingBall')

                %%%Load and Sort directions

                if a == 3

                    stimDir =convertStringsToChars(string(basicPathA{a})+"\"+string(expA{a})+"\Insertion"+(in+1));
                else

                    stimDir =convertStringsToChars(string(basicPathA{a})+"\"+string(expA{a})+"\Insertion"+in);
                end

                file = dir (stimDir);
                filenames = {file.name};
                ballFiles = filenames(contains(filenames,"linearlyMovingBall"));
                directions = [];
                offsets = [];

                j =1;
                if size(ballFiles) ~= [0 0]

                    for i = ballFiles
                        ball= load(stimDir+"\"+string(i));

                        directions = [directions cell2mat(ball.VSMetaData.allPropVal(17))];

                        offsets = [offsets cell2mat(ball.VSMetaData.allPropVal(18))];

                        direcNames = unique(directions);

                        stimDurStats = cell2mat(ball.VSMetaData.allPropVal(38))*1000;
                        interStimStats = cell2mat(ball.VSMetaData.allPropVal(28))*1000;

                        j = j+1;
                    end
                    disp('Visual stats extracted!')
                else
                    disp('Directory does not exist!');
                end

                             
                A = [stimOn' directions' offsets'];

                C = sortrows(A,[2 3]);

                %Sort directions:

                directimesSorted = C(:,1)';

                [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((stimDur+preBase)/bin));

                %%%% Convolute matrix:

                %%%1. Convolute in the 3rd dimension (trials)

                [MrC]=ConvBurstMatrix(Mr,fspecial('gaussian',[1 5],3),'same');


                [nT,nN,nB] = size(MrC);

                %%%2.Convolute in the 1st dimension (bin) and reduce
                %%%element number to number of offsets
                
                trialDivision = nT/(length(unique(offsets))*length(unique(directions)));

                meanDirOff = zeros(nT/(trialDivision), nN, nB);

                j = 1;

                for i = 1:trialDivision:nT


                    meani = mean(MrC(i:min(i+trialDivision-1, end),:,:),1);

                    meanDirOff(j,:,:) = meani;

                    j = j+1;

                end


                %%% end of convolution



                %%%%Create window of 2 offsets by 10 bins (500 ms) to scan
                %%matrix and select highest mean druing scan.

                % Define the size of the window
                if a == 1 %PV67
                    window_size = [4*trialDivision, 2]; %bin = 50 ms, 4 offsets
                else
                    window_size = [trialDivision, 2];
                end


                % Initialize the maximum mean value and its position %No
                % need or concolutions. 

                max_position = zeros(nN,2);
                max_mean_value = zeros(1,nN);

                for u =1:nN
                    % Slide the window over the matrix
                    %unit matrix
                    max_mean_value(u) = -Inf;
                    for d = 1:length(unique(offsets))*trialDivision:nT
                        uM = squeeze(MrC(d:d+length(unique(offsets))*trialDivision-1,u,:));
                        for i = 1:size(uM, 1) - window_size(1) + 1
                            for j = 1:size(uM, 2) - window_size(2) + 1
                                % Extract the sub-matrix
                                sub_matrix = uM(i:min(i+window_size(1)-1, end), j:min(j+window_size(2)-1,end));

                                % Compute the mean value
                                mean_value = mean(sub_matrix(:));

                                % Update the maximum mean value and its position
                                if mean_value > max_mean_value(u)
                                    max_mean_value(u) = mean_value;
                                    max_position(u,:) = [i+d-1, j];

                                end


                            end
                        end
                    end
                end

                spkRateR = max_mean_value;

                epsilon = 0.01;

                denom = mad(MbC,0)+epsilon; %mean(Mb,0)+epsilon; %

                mSpk = mean(spkRateR);

                Zscore = (spkRateR - (spkRateBM + spkRateR)/2)./denom;

                RespU = GoodU_or(Zscore>2);


                sDirections = sort(directions);

                tunning = sDirections(max_position(:,1));



                %MinusRB = (mean(spkRateR) -(mean(spkRateB)))./(mean(spkRateR) + (mean(spkRateB)));

                MinusRB = spkRateR - spkRateBM;


                %                 if a==2 & in ==1
                %                     2+2
                %
                %                 end

                if (a == 1 && in ==1) || (a == 2 && in ==1) || (a == 3 && in ==1)
                    plotexamplesMB =0;

                    if plotexamplesMB == 1
                        %                     % Create direction arrows:
                        %
                        %                     % Define the x and y coordinates for the arrows
                        %                     x = zeros(1,direcN);
                        %                     y = categ/direcN:categ/direcN:categ;
                        %                     theta = unique(directions);
                        %
                        %                     % Calculate the u and v components of the arrows
                        %                     u = cos(theta);
                        %                     v = sin(theta);
                        %
                        %                     % Normalize the arrow vectors to make them all the same length
                        %                     arrow_length = 2;
                        %                     u = u ./ sqrt(u.^2 + v.^2) * arrow_length;
                        %                     v = v ./ sqrt(u.^2 + v.^2) * arrow_length;
                        %
                        %                     % Adjust the starting points of the arrows
                        %                     x = x - u/2;
                        %                     y = y - v/2;
                        %
                        %                     figd = figure;
                        %                     h = quiver(x, y, u, v, 'AutoScale', 'off');
                        %                     set(h, 'Color', 'k', 'LineWidth', 1);% Set the aspect ratio% Set the color to black and the line width to 2
                        %                     set(gca, 'DataAspectRatio', [1 1 1]);
                        %                     set(gca, 'XTick', [], 'YTick', [], 'XColor', 'none', 'YColor', 'none'); % Remove ticks and axis borders
                        %                     axis([-1 4 0 40]);
                        %                     set(figd, 'Color', 'white');
                        %                     print(figd, sprintf('%s-In.%d-MovBall-dir.png',expA{a},in),'-dpng');
                        %
                        %                     close all


                        for u = eNeuron
                            [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round((stimDur+preBase*2)/bin));

                            [nT,nN,nB] = size(Mr);

                            fig = figure;
                            categ = nT/trialDivision;
                            offsetN = length(unique(offsets));
                            direcN = length(unique(directions));
                            imagesc(squeeze(Mr(:,u,:)));colormap(flipud(gray(64)));
                            xline(preBase/bin,'k', LineWidth=1.5)
                            xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5)
                            ylabel('Trials');xlabel('Time (ms)');
                            m = max_mean_value(u) - spkRateBM(u);
                            title(sprintf('U.%d-R.%.3f-B.%.3f-S.%.3f',u,max_mean_value(u),spkRateBM(u),Zscore(u)));
                            xticks([0.5 (preBase/bin):10:nB])
                            xticklabels([-preBase 0:10*bin:nB*bin])
                            v = nT/direcN:nT/direcN:nT-1;
                            yline(v+0.5,'r--', LineWidth=3);
%                             hcb = colorbar();
%                             title(hcb,'Spikes/sec');
                            caxis([0 max(0.2,max(max_mean_value(u)))])
                            hold on
                            rectangle('Position', [max_position(u,2), max_position(u,1), window_size(2), window_size(1)], 'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');
                            hold off
                            prettify_plot
                            %fig.Position = [-755   440   520   438];
                            print(fig, sprintf('%s-In.%d-MovBall-dir-U%d-W%d-%dW.png',expA{a},in,u,window_size(1),window_size(2)),'-dpng');
                            close

                        end

                    end
                end



            end

            %spkRateR = sum(Mr,3)./(stimDur/1000);

            

            MinusRBs{1,s} = MinusRB;
            MinusRBs{2,s} =stimName;% MinusRB;

            Zscores{1,s} = Zscore;
            Zscores{2,s} = stimName;

            tune{1,s} = tunning;
            tune{2,s} = stimName;


        end

        MinusRBi{in} = MinusRBs;
        Zscorei{in} = Zscores;
        Tunei{in} = tune;


        itN = itN+1;

        %         respIns = cell2mat(MinusRBs(1,:));
        %
        %         q1 = quantile(respIns, 0.1);
        %
        %         q3 = quantile(respIns, 0.9);
        %
        %         MinusRBi{2,in} = [q1 q3];

    end

    MinusRBt{1,a} = MinusRBi;

    MinusRBt{2,a} = Zscorei;

    MinusRBt{3,a} = Tunei;

    MinusRBt{4,a} = expA{a};


end

%% Responsive neurons to Mov ball vs response in rect grid. Change in spike rate


valsMBzs = [];
valsMBsr = [];


valsRGzs= [];
valsRGsr= [];


valsFFzs= [];
valsFFsr= [];

valsSDzs = [];
valsSDsr = [];

valsRG = [];
valsMB = [];
valsSD = [];
valsFF = [];


for a = 1:length(expA)

    for in =  1:length(animalA{a})

        for s = 1:3


            stimName = StimOrder{stimTTL{a}(2,s)};


            if string(stimName) == "linearlyMovingBall"

                valsZS = MinusRBt{2,a}{1,in}{1,s};

                valsSR = MinusRBt{1,a}{1,in}{1,s};

                indMB = find(valsZS>2);

                valsMBzs = [valsMBzs valsZS(indMB)];

                valsMBsr = [valsMBsr valsSR(indMB)];

                valsMB = [valsMB valsZS];

                

            end

        end

        for s = 1:3

            stimName = StimOrder{stimTTL{a}(2,s)};

            valsSZ = MinusRBt{2,a}{1,in}{1,s};
            valsSR = MinusRBt{1,a}{1,in}{1,s};



            if string(stimName) == "rectGrid"

                valsRG = [valsRG valsSZ];

                valsRGzs = [valsRGzs valsSZ(indMB)];

                valsRGsr =  [valsRGsr valsSZ(indMB)];


            end

            if string(stimName) == "fullFieldFlash"

                valsFF = [valsFF valsSZ];

                valsFFzs = [valsFFzs valsSZ(indMB)];

                valsFFsr = [valsFFsr valsSZ(indMB)];

            end

            if string(stimName) == "StaticDriftingGrating"

                valsSD = [valsSD valsSZ];

                valsSDzs  = [valsSDzs valsSZ(indMB)];

                valsSDsr = [valsSDsr valsSZ(indMB)];

            end

        end

    end

end

SZs =[valsMBzs,valsRGzs,valsFFzs,valsSDzs];

SRs = [valsMBsr, valsRGsr,valsFFsr,valsSDsr];

MBvxRG = [valsMB valsRG];

nRmb = length(valsMB(valsMB>2));

nRrg = length(valsRG(valsRG>2));

grp = [zeros(1,length(valsMBzs)),ones(1,length(valsRGzs)), ones(1,length(valsFFzs))+1, ones(1,length(valsSDzs))+2];

grp2 = [zeros(1,length(valsMB)),ones(1,length(valsRG))];

grp3 = [zeros(1,length(valsMBzs)),ones(1,length(valsRGzs))];


sFF = valsFF(valsFF>2);
sSD = valsSD(valsSD>2);


%%Violin plot
figure
violinplot(SZs,grp,'QuartileStyle','shadow','ShowMedian', true,'DataStyle', 'histogram');
ylabel('Z-score')
ylim([0 10])
prettify_plot


f = plotSpread(datac);
set(findall(1,'type','line'),'markerSize',16);
ylabel('Z-score')
ylim([0 10])
prettify_plot
set(f{1},'color',[0 0.3 0.6]);
set(gcf, 'Color', 'white');




%

figure
histogram(valsMB,'BinEdges',[0:0.5:10])
xlim([0 10])
hold on
histogram(valsRG,'BinEdges',[0:0.5:10])
xlim([1.5 10])
ylim([0 250])
xline(2,'r',LineWidth=1.5)
xlabel('Z-score')
ylabel('Unit count')
legend({'Moving ball (697 r.u.)','Grid of rectangles (133 r.u.)'})
prettify_plot




figure()

h = boxplot(MBvxRG, grp2, 'Notch', 'on','Labels',{'LMB','RG'}, 'ExtremeMode','clip',...
    'Colors',[0.5 0.7 0.7]);
set(h, 'LineWidth', 2); 
grid on
hold on
%f = plotSpread(datac);
set(findall(1,'type','line'),'markerSize',16);
ylabel('Z-score')
ylim([0 10])
prettify_plot
% set(f{1},'color',[0 0.3 0.6]);
% set(gcf, 'Color', 'white');


valstimN = (valStim - mean(valStim))./std(valStim);

maxVSN = max(valstimN);

minVSN = min(valstimN);

MaxVal = ceil(MaxVal*10)/10;
MinVal = floor(MinVal*10)/10;

%Mvalues = round(normalize(linspace(MinVal,MaxVal,50),'range',[10 200]));

%%
meanValsIns = {MrFFF,MrLMB,MrRG,MrSDG};

names = {'fullFieldFlash','linearlyMovingBall','rectGrid','StaticDriftingGrating'};

figure;
t = tiledlayout(2,2,"TileSpacing","compact");

for v = 1:length(meanValsIns)

    nexttile

    non_zero_rows = any(meanValsIns{v}, 2); % Logical index for non-zero rows
    non_zero_columns = any(meanValsIns{v}, 1); % Logical index for non-zero columns


    fmatrix =meanValsIns{v}(non_zero_rows,non_zero_columns);

    %matrix = reshape(matrix,length(matrix)/size(meanValsIns{v}(meanValsIns{v}(1,:)~=0), 2),size(meanValsIns{v}(meanValsIns{v}(1,:)~=0),2));
    
    y = mean(fmatrix);

    std_errs = std(fmatrix, 0, 1)/sqrt(size(fmatrix, 1)); 

    x = (1:length(y))*bin;

    ysmooth = smooth(x,y);

%     % Generate a finer x-axis for smooth plotting
%     x_interp = linspace(min(x), max(x), 1000);
% 
%     % Perform cubic spline interpolation
%     y_interp = interp1(x, y, x_interp, 'spline');

    plot(x,ysmooth,'LineWidth',2)
    hold on
    plot(x,y,'.','Color','b','MarkerSize',20)
    title(names{v})
    xline([preBase2 max(x)-preBase2],'LineWidth',2)
    %errorbar(x, y, std_errs, 'o', 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b')
    xlim([0,max(x)])
    xticks(0:preBase2:max(x))
    xticklabels(-preBase2:preBase2:max(x)-preBase2)

    if length(names{v}) == length('StaticDriftingGrating')

        xline([4*bin+preBase2 12*bin+preBase2],'Color','r','LineWidth',2)


    elseif length(names{v}) == length('fullFieldFlash')

         xline([2*bin+preBase2 8*bin+preBase2],'Color','r','LineWidth',2)


    elseif length(names{v}) == length('rectGrid')

         xline([2*bin+preBase2 8*bin+preBase2],'Color','r','LineWidth',2)


    elseif length(names{v}) == length('linearlyMovingBall')

         xline([4*bin+preBase2 12*bin+preBase2],'Color','r','LineWidth',2)
    end



    hold off

end
%linkaxes(findobj(gcf, 'Type', 'Axes'), 'y');

ylabel(t, 'Mean response (spikes/sec)');

xlabel(t,'Time (ms)');


prettify_plot

cd('\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\')
save("ResponsesAnimals","MinusRBt")

