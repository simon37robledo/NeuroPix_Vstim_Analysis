%function MovBallRasters(NP, ball, in, ttl_index)
%% Phy commmand: phy template-gui params.py
%Order of batch= FFF, movBall, rectGrid, RectNoise, Gratings.

%% PROCESSING MOVING BALL STIMULUS

ins = 1; %detect how many ins
basic_path = '\\132.66.45.127\data\Large_scale_mapping_NP\\Immobilized_exp\PV67';
exp = 'PV67_experiment_5_7_23';

animal = "PV67";
ttl_index =2;
in_ind =1; %initialize insertion index. 

animalSum = cell(1,length(ins));

allDepth =  cell(1,length(ins));

allXdist =  cell(1,length(ins));

DirSpkR = cell(1,length(ins));

DirSpkB = cell(1,length(ins));

DirSpkZB = cell(1,length(ins));

DirSpkZA = cell(1,length(ins));

GoodUnits = cell(1,length(ins));

minusSpK = cell(1,length(ins));


rasters = 1;
heatmaps = 0;
heatmap_fig = 0;
NumBatch = 0; 

%InsAngles = [88 88]; %PV102
% 
%InsDepths = [3909 3931];%PV102

InsDepths = [3913.9 3904.34 3525.7 3900.1 3914.34 3739.3 3906.8];

% InsAngles = [88 88 88 88];%PV67
% 
% InsDepths = [3956 3907 4000 3934]; %PV67

%%

for in = ins %Start iterating across insertions



    path = convertStringsToChars(basic_path+string(filesep)+exp+string(filesep)+"Insertion"+in+string(filesep)+"catgt_"+exp+"_"+in+"_g0");

    NP = NPAPRecording(path);

    [Ttrigger,chNumberT]=NP.getTrigger();

    Sstart = [Ttrigger{1}];
    Send = [Ttrigger{2}];
    if animal == "PV27" & length(Send)>3
         Send = Send(2:end); %PV27 weirdness, it has one more ttl off trigger at the start
    end
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

    %--------------------------------------------------------------------------------------------------------------------------------%
    % Load spike sorting and Phy results

    % Convert into format to generate raster plots:
    p = NP.convertPhySorting2tIc(NP.recordingDir);

    %Select good units
    label = string(p.label');

    goodU = p.ic(:,label == 'good');
    amp = p.neuronAmp(label=='good');

    disp("Phy extracted!")

    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');



    cluster_group = readtable(string(NP.recordingDir) + "\cluster_group.tsv",  "FileType","text",'Delimiter', '\t');


    %--------------------------------------------------------------------------------------------------------------------------------%
    % Build raster plots per unit
    bin=20;
    win=stimDurStats+interStimStats;
    start = interStimStats/2;


    %image saving path:
    pathF =  convertStringsToChars(basic_path+"\"+exp+"\Insertion"+in+"\Figures\MovBall");
    cd(pathF)

    A = [stimOn' directions' offsets'];

    C = sortrows(A,[2 3]);

    %Sort directions:

    directimesSorted = C(:,1)';

    unitCh = zeros(length(goodU(1,:)),2);

    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");


    %--------------------------------------------------------------------------------------------------------------------------------%
    % Build raster plots per unit
    GoodUnits{in_ind} = GoodU_or;

    if rasters ==1 
        
     
        for u = 1:length(GoodU_or)
            % Select stimulus
            %u = 210; %WEST BOI

            %u = 121; 

            unitCh(u,:) =cluster_info.ch(cluster_info.cluster_id == GoodU_or(u));

            direcNames = unique(directions);

            offsetNames = unique(offsets);

            stim = directimesSorted - start;

            %j = NP.getData(130,500,1000);

            title_stim = sprintf('Insertion-%d-Directions-Unit(%d)-%d-channel-#%d', in, GoodU_or(u),u, cluster_info.ch(cluster_info.cluster_id == GoodU_or(u)));
                
            bin = 20;
            [M]=BuildBurstMatrix(goodU,round(p.t/bin),round(stim/bin),round(win/bin));

            sum(M(:,u,:),'all')

            [nTrials,nNeurons,nTimes]=size(M);

            categ = (length(direcNames)*length(offsetNames))/2;

            categ_name = 1:categ;
                
            if offsets ==1
            
                Z = [];

                direcNamesMod = string(direcNames);

                offsetNamesMod = string(offsetNames);

                direcNamesMod(direcNamesMod == "0") = "North";

                direcNamesMod(direcNamesMod == "1.5708") = "West";

                direcNamesMod(direcNamesMod == "3.1416") = "South";

                direcNamesMod(direcNamesMod == "4.7124") = "East";

                for i = 1:length(direcNames)
                    for j = 1:length(offsetNames)
                        Z = [Z, strcat(direcNamesMod(i), " & ", string(offsetNames(j)))];
                    end
                end
                direcOffsets = length(Z)/length(direcNames);

                Z1 = Z([1:direcOffsets direcOffsets*2+1:direcOffsets*3]);
                Z2 = Z([direcOffsets+1:direcOffsets*2 direcOffsets*3+1:end]);


                nTrials = nTrials/2;

                direcTrials = length(M)/length(direcNames);

                h = figure;
                imagesc((1:nTimes)*bin,1:2:2*nTrials,squeeze(M([1:direcTrials direcTrials*2+1:direcTrials*3],u,:)));colormap(flipud(gray(64))); %Vertical directions
                ylabel('Trials');xlabel('Time (ms)');
                title(sprintf('%s -%s trials per category-Vertical mov.',title_stim, string(nTrials/categ)));
                %cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');clim([0 1]);
                xline(start, '-g', 'Start','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
                xline(start+stimDur, '-b', 'End', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
                %xline(2000, '-r', 'Inter-trial', 'LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom');
                v = 1:(2*nTrials)/categ:2*nTrials;
                yline(v(length(v)/2+1), 'r', LineWidth=2);
                yline(v, '--k', Z1, 'LabelVerticalAlignment','bottom', LineWidth=1);
                yline(v(length(v)/2+1), 'r', LineWidth=2);
                x = 1:2*nTrials/categ:2*nTrials; %to add first ticks (even);
                y = nTrials/categ:2*nTrials/categ:2*nTrials; %add second ticks (uneven);
                yticklabels(repmat([1  round((nTrials/categ)/2)],1,categ)); yticks(sort([x y]));
                caxis([0 1]);
                set(gcf, 'units','normalized','outerposition',[0 0 0.6 0.6]);
                orient(gcf,'landscape');
                pause(0.5)
                print(h, sprintf('%s-vertical',title_stim),'-dpng');
                %print(h, sprintf('%s-vertical',title_stim),'-dpdf','-vector');
                prettify_plot

                clearvars h
                %close all


                h = figure;
                imagesc((1:nTimes)*bin,1:2:2*nTrials,squeeze(M([direcTrials+1:direcTrials*2 direcTrials*3+1:end],u,:)));colormap(flipud(gray(64))); %Vertical directions
                ylabel('Trials');xlabel('Time (ms)');
                title(sprintf('%s -%s trials per category-Horizontal mov.',title_stim, string(nTrials/categ)));
                %cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');clim([0 1]);
                xline(start, '-g', 'Start','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
                xline(start+stimDur, '-b', 'End', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
                %xline(2000, '-r', 'Inter-trial', 'LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom');
                v = 1:(2*nTrials)/categ:2*nTrials;
                yline(v(length(v)/2+1), 'r', LineWidth=2);
                yline(v, '--k', Z2, 'LabelVerticalAlignment','bottom', LineWidth=1);
                x = 1:2*nTrials/categ:2*nTrials; %to add first ticks (even);
                y = nTrials/categ:2*nTrials/categ:2*nTrials; %add second ticks (uneven);
                yticklabels(repmat([1  round((nTrials/categ)/2)],1,categ)); yticks(sort([x y]));
                caxis([0 1]);
                set(gcf, 'units','normalized','outerposition',[0 0 0.6 0.6]);
                orient(gcf,'landscape');
                pause(0.5)
                print(h, sprintf('%s-horizontal',title_stim),'-dpng');
                %print(h, sprintf('%s-horizontal',title_stim),'-dpdf','-vector');
                prettify_plot

                clearvars h
                close all
            else
                h = figure;
                categ = length(direcNames);
                imagesc((1:nTimes)*bin,1:2:2*nTrials,squeeze(M(:,u,:)));colormap(flipud(gray(64))); %Vertical directions
                ylabel('Trials');xlabel('Time (ms)');
                title(sprintf('%s -%s trials per category.',title_stim, string(nTrials/categ)));
                %cb=colorbar;set(cb,'Position',[0.91 0.75 0.013 0.17]);ylabel(cb,'Spk. Rate');clim([0 1]);
                xline(start, '-g', 'Start','LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom', LineWidth=1.5);
                xline(start+stimDur, '-b', 'End', 'LabelHorizontalAlignment', 'right', LabelVerticalAlignment='bottom', LineWidth=1.5);
                %xline(2000, '-r', 'Inter-trial', 'LabelHorizontalAlignment', 'right',LabelVerticalAlignment='bottom');
                v = 1:(2*nTrials)/categ:2*nTrials;
                yline(v, '--k');
                x = 1:2*nTrials/categ:2*nTrials; %to add first ticks (even);
                y = nTrials/categ:2*nTrials/categ:2*nTrials; %add second ticks (uneven);
                yticklabels(repmat([1  round((nTrials/categ)/2)],1,categ)); yticks(sort([x y]));
                caxis([0 1]);
                % Get the current axes
                ax = gca;

                % Add the text
                for t = 1:length(v)
                    text(ax.XLim(2) + 50, v(t), string(direcNames(t)), 'VerticalAlignment', 'top');
                end

                % Adjust the x-axis limit to make room for the text
                ax.XLim = [ax.XLim(1) ax.XLim(2) + 0.2];
                set(gcf, 'units','normalized','outerposition',[0 0 0.6 0.6]);
                orient(gcf,'landscape');
                pause(0.5)
                prettify_plot
                print(h, sprintf('%s-horizontal',title_stim),'-dpng');
                print(h, sprintf('%s-horizontal',title_stim),'-dpdf', '-bestfit','-vector');


                clearvars h
                close all

            end

            writematrix(unitCh,sprintf('In-%d-units-chan.xls',in))
        end

    %--------------------------------------------------------------------------------------------------------------------------------%
    % Position

    %Spike index and times:

    spkI = readNPY(string(NP.recordingDir) + "\spike_clusters.npy");

    spkt = readNPY(string(NP.recordingDir) + "\spike_times.npy")/(NP.samplingFrequency/1000);

    %Good units

    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

    %Insertion info:
    AngleInser = InsAngles(in_ind);
    InserDepth = InsDepths(in_ind);
    ShankDist = zeros(1, length(GoodU_or)); %distance of unit from surface along the shank
    verticalDepth = zeros(1, length(GoodU_or)); %depth of unit along vertical axis
    XDist = zeros(1, length(GoodU_or)); %X distance of unit from insertion


    %%%

    %Get X and Y positions

    X = squeeze(cell2mat(ball.VSMetaData.allPropVal(22)));

    Y = squeeze(cell2mat(ball.VSMetaData.allPropVal(23)));

    %times values of frames per trial

    times = linspace(0,round(stimDur),length(X));

    sX = size(X);

    %number of trials

    trials = cell2mat(ball.VSMetaData.allPropVal(29));

    %%%



    Xp =  permute(X, [3 1 2]); %Change order of dimensions so the reshaped vector has all frames-Ofsset1-Dir1, then frames-offset2-dir1.... then frames-offset1-dir2
    %Xtest = reshape(Xp, 1, []);
    Xv =  reshape(Xp, 1, []);


    n = numel(Xv);

    direcDiv = length(Xv)/length(unique(directions));
    Xvd = mat2cell(Xv,1,diff([0:direcDiv:n-1,n]));

    Yp = permute(Y, [3 1 2]);
    Yv = reshape(Yp, 1, []);

    Yvd = mat2cell(Yv,1,diff([0:direcDiv:n-1,n]));

    %Count istances of coordinate pairs.

    hozD = unique(round(Y(:,2,:))); %horizontal direction (east)

    vertD = unique(round(X(:,1,:))); %vertical direction (north)

    binWidth = diff(hozD(1:2));

    EdgesY = [hozD(1)-binWidth/2 ; hozD+binWidth/2]; %Create edges Y

    EdgesX = [vertD(1)-binWidth/2 ; vertD+binWidth/2]; %create edges X

    counts = cell(1,length(Xvd));

    for i = 1:length(Xvd)


        counts{1,i} = histcounts2(cell2mat(Xvd(1,i)), cell2mat(Yvd(1,i)), EdgesX, EdgesY, 'Normalization','count'); % Count instances if the ball in each coordinate

    end
    %imagesc(EdgesX, EdgesY, counts); %plot it




    %Xpos = unique(round(X/a)); %unique X positions

    %Ypos= unique(round(Y/a)); %unique Y positions


    rx = [];
    ry = [];
    t1 = [];

    %Create

    vec = 1:length(times):length(Xv)+1;

    for i =1:length(vec)-1

        rx = [rx repmat(Xv(vec(i):vec(i+1)-1), 1, trials)]; %15 = trials.

        ry = [ry repmat(Yv(vec(i):vec(i+1)-1), 1 ,trials)];

        t1 = [t1 repmat(times, 1, trials)];


    end



    test2 = Xv(length(Xv)-length(times)*2:length(Xv)-length(times));

    test1 = rx(length(rx)-length(times)*16:length(rx)-length(times)*15);

    j =1;


    for i= 1:length(times):length(t1) %iterate trough every trial and assign times


        t1(i:i+length(times)-1) = t1(i:i+length(times)-1)+directimesSorted(j);


        j = j+1;

    end


    sumNeurons = zeros(4,length(GoodU_or));

    frameDur = stimDur/length(times);

    baseLine = 500;
    win = stimDur+baseLine;


    [M]=BuildBurstMatrix(goodU,round(p.t/frameDur),round((directimesSorted-baseLine)/frameDur),round(win/frameDur));

    stdBL = std(M(:,:,1:round(baseLine/frameDur)),0,3);

    meanBL = mean(M(:,:,1:round(baseLine/frameDur)),3);

    stdBLperN = std(meanBL);

    meanBLperN = mean(meanBL);


    %replace zero values

    stdBL(stdBL==0) = mean(stdBLperN);

    stdBLperN(stdBLperN==0) = mean(stdBLperN);


    %MnT=bsxfun(@rdivide,bsxfun(@minus,M,meanBL),stdBL); %normalization per triak

    %MnN =bsxfun(@rdivide,bsxfun(@minus,M,meanBLperN),stdBLperN); %normalization per neuron


    %allM = cell(4,length(GoodU_or));

    allM1 = zeros(length(GoodU_or),7,7);
    allM2 = zeros(length(GoodU_or),7,7);
    allM3 = zeros(length(GoodU_or),7,7);
    allM4 = zeros(length(GoodU_or),7,7);

    %ITERATE ACROS UNITS

    for u = 1:length(GoodU_or)

        %u = 181;

        %Main unit channel position

        %change channel 0 to 1

        cluster_info.ch(cluster_info.ch ==0) = 1;

        ch = cluster_info.ch(cluster_info.cluster_id == GoodU_or(u));

        ShankDist(u) = InserDepth - (NP.chLayoutPositions(2,ch));

        verticalDepth(u) = sin(deg2rad(AngleInser))*ShankDist(u); %depth of unit along vertical axis

        XDist(u) = cos(deg2rad(AngleInser))*ShankDist(u); %X distance of unit from insertion



        title_stim = sprintf('Insertion-%d-Unit-%d-channel-#%d', in, GoodU_or(u), ch);


        Unit_spks = spkt(spkI == GoodU_or(u));


        if heatmaps ==1

            %     MuT = squeeze(MnT(:,u,round(baseLine/frameDur)+1:end));
            %
            %     MuT = reshape(MuT,1,[]);
            %
            %
            %     MuN = squeeze(MnN(:,u,round(baseLine/frameDur)+1:end));
            %
            %
            %     MuN = reshape(MuN ,1,[]);
            %
            %     sC = MuN;


            sC =[];

            for i= 1:length(t1) %iterate every time and asign number of spikes


                sC = [sC sum(Unit_spks > t1(i) & Unit_spks < t1(i)+(t1(2)-t1(1)))/(frameDur)];


            end

            BL = zeros(1,length(t1));



            for i=1:sX(3):length(t1) %Iterate to get baseline (1000 ms before start)

                BL(i:i+sX(3)-1) = sum((Unit_spks > t1(i)-interStimStats/4 & Unit_spks < t1(i)))/(baseLine);

                %sC(i:i+sX(3)-1) = (sC(i:i+sX(3)-1) - BL(i:i+sX(3)-1))/std(BL(i:i+sX(3)));

            end

            B = [t1;rx;ry;sC;BL]; %Create matrix

            %separate matrix into 2
            B1 = B(:,1:length(B)/4);  %North

            B2 = B(:,length(B)/4+1:length(B)/2);%West

            B3 = B(:,length(B)/2+1:(length(B)*3/4)); %South

            B4 = B(:,(length(B)*3/4)+1:end); %East

            M1 = zeros(length(hozD),length(vertD));
            M1b = zeros(length(hozD),length(vertD));

            M2 = zeros(length(hozD),length(vertD));
            M2b = zeros(length(hozD),length(vertD));

            M3 = zeros(length(hozD),length(vertD));
            M3b = zeros(length(hozD),length(vertD));

            M4 = zeros(length(hozD),length(vertD));
            M4b = zeros(length(hozD),length(vertD));

            %  [countsHalfNS] = histcounts2(B1(2,:), B1(3,:), EdgesX, EdgesY, 'Normalization','count');
            %
            %  [countsHalfEW] = histcounts2(B2(2,:), B2(3,:), EdgesX, EdgesY, 'Normalization','count');

            Rms = 100; %response period in ms


            frameMs = stimDur/length(times);

            framResp = round(Rms/frameMs);

            %%Normalize std:



            for x =1:length(hozD)

                for y = 1:length(vertD) %How to select the spikes that are in a position between the edges?

                    %%% North

                    cond1 = (B1(2,:) > vertD(x) - binWidth/2) & (B1(2,:) < vertD(x) + binWidth/2);

                    cond2 = (B1(3,:) >  hozD(y) - binWidth/2) & (B1(3,:) <  hozD(y) + binWidth/2);

                    cols1 = cond1 & cond2;

                    %M1(x,y) = mean((B1(4,cols1)-mean(B1(5,cols1)))/std(B1(5,cols1)));
                    M1(x,y) = mean(B1(4,cols1));
                    M1b(x,y) = mean(B1(5,cols1));


                    %%% West

                    cond3 = B2(2,:) > vertD(x) - binWidth/2 & B2(2,:) < vertD(x) + binWidth/2;

                    cond4 = B2(3,:) >  hozD(y) - binWidth/2 & B2(3,:) <  hozD(y) + binWidth/2;

                    cols2 = cond3 & cond4;
                    %M2(x,y) = mean((B2(4,cols2)-mean(B2(5,cols2)))/std(B2(5,cols2)));
                    M2(x,y) = mean(B2(4,cols2));
                    M2b(x,y) = mean(B2(5,cols2));


                    %%% South

                    cond5 = B3(2,:) > vertD(x) - binWidth/2 & B3(2,:) < vertD(x) + binWidth/2;

                    cond6 = B3(3,:) >  hozD(y) - binWidth/2 & B3(3,:) <  hozD(y) + binWidth/2;

                    cols3 = cond5 & cond6;
                    %M3(x,y) = mean((B3(4,cols3)-mean(B3(5,cols3)))/std(B3(5,cols3)));
                    M3(x,y) = mean(B3(4,cols3)); %mean instead of sum???
                    M3b(x,y) = mean(B3(5,cols3));


                    %%% East

                    cond7 = B4(2,:) > vertD(x) - binWidth/2 & B4(2,:) < vertD(x) + binWidth/2;

                    cond8 = B4(3,:) >  hozD(y) - binWidth/2 & B4(3,:) <  hozD(y) + binWidth/2;

                    cols4 = cond7 & cond8;
                    %M4(x,y) = mean((B4(4,cols4)-mean(B4(5,cols4)))/std(B4(5,cols4)));
                    M4(x,y) = mean(B4(4,cols4));
                    M4b(x,y) = mean(B4(5,cols4));

                end
            end


            %     %%Normalize
            %     nM1 = rdivide(M1,countsHalfNS);
            %
            %     nM2 = rdivide(M2,countsHalfEW);
            %
            %     nM3 = rdivide(M3,countsHalfNS);
            %
            %     nM4 = rdivide(M4,countsHalfEW);


            %%Normalize by dividing by the max value.
            %maxVal = unique(max(cat(2, M1,M2,M3,M4),[],'all'));



            %     nM1 = (M1/maxVal)';
            %
            %     nM2 = (M2/maxVal)';
            %
            %     nM3 = (M3/maxVal)';
            %
            %     nM4 = (M4/maxVal)';

            concatMb = cat(2, M1b,M2b,M3b,M4b);

            stdBLrect = std(concatMb,0,'all',"omitnan");

            stdBLrect(stdBLrect==0) = stdBLperN(u);

%             [MbAll] = BuildBurstMatrix(goodU,round(p.t),round(stimOn-stimInter),round(stimInter));
% 
%             spkRateBall = sum(squeeze(Mb(:,u,:)),2)/(stimInter/1000);
% 
%             epsilon = 0.1;
% 
%             stdBLrect = std(spkRateBall)+epsilon;

            %
            nM1 = ((M1-mean(concatMb,'all',"omitnan"))/stdBLrect);

            nM2 = ((M2-mean(concatMb,'all',"omitnan"))/stdBLrect);

            nM3 = ((M3-mean(concatMb,'all',"omitnan"))/stdBLrect);

            nM4 = ((M4-mean(concatMb,'all',"omitnan"))/stdBLrect);

            nM1 = nM1(2:end-1,2:end-1);
            nM2 = nM2(2:end-1,2:end-1);
            nM3 = nM3(2:end-1,2:end-1);
            nM4 = nM4(2:end-1,2:end-1);

            %     nM1 = (M1)';
            %
            %     nM1(isinf(nM1)|isnan(nM1)) = 0;
            %
            %     nM2 = (M2)';
            %
            %     nM2(isinf(nM2)|isnan(nM2)) = 0;
            %
            %     nM3 = (M3)';
            %
            %     nM3(isinf(nM3)|isnan(nM3)) = 0;
            %
            %     nM4 = (M4)';
            %
            %     nM4(isinf(nM4)|isnan(nM4)) = 0;

            maxVal = unique(max(cat(2, nM1,nM2,nM3,nM4),[],'all'));
            minVal =  unique(min(cat(2, nM1,nM2,nM3,nM4),[],'all'));
            concatMn = cat(2, nM1,nM2,nM3,nM4);

            %%Plot

            if heatmap_fig

                t = tiledlayout(2, 2,'TileSpacing','tight');

                LB=flipud(lbmap(256,'BrownBlue'));

                axisL = 1:length(nM1);

                nexttile
                imAlpha=ones(size(nM1));
                imAlpha(isnan(nM1))=0;
                imagesc(nM1,'AlphaData',imAlpha);
                %plot it
                title('North');
                colormap(LB);
                caxis([-max(abs(concatMn(:)))-0.1 max(abs(concatMn(:)))]);
                ax = gca;
                axis equal
                axis tight
                %set(gca,'color',0*[1 1 1]);
                set(gca,'XTick',[],'YTick',[]);

                nexttile
                imAlpha=ones(size(nM3));
                imAlpha(isnan(nM3))=0;
                imagesc(nM3,'AlphaData',imAlpha); %plot it
                title('South');
                colormap(LB);
                caxis([-max(abs(concatMn(:)))-0.1 max(abs(concatMn(:)))]);
                ax = gca;
                axis equal
                axis tight
                %set(gca,'color',0*[1 1 1]);
                set(gca,'XTick',[],'YTick',[]);
                % ax.CLim = [0  maxVal];
                cb = colorbar();
                cb.Ticks = [-30, -15, 0, 15, 30];


                nexttile
                imAlpha=ones(size(nM4));
                imAlpha(isnan(nM4))=0;
                imagesc(nM4,'AlphaData',imAlpha); %plot it
                title('East');
                colormap(LB);
                caxis([-max(abs(concatMn(:)))-0.1 max(abs(concatMn(:)))]);
                ax = gca;
                axis equal
                axis tight
                %set(gca,'color',0*[1 1 1]);
                set(gca,'XTick',[],'YTick',[]);
                %ax.CLim = [0  maxVal];
                %colorbar('gray');


                nexttile
                imAlpha=ones(size(nM2));
                imAlpha(isnan(nM2))=0;
                imagesc(nM2,'AlphaData',imAlpha); %plot it
                title('West');
                colormap(LB);
                caxis([-max(abs(concatMn(:)))-0.1 max(abs(concatMn(:)))]);
                ax = gca;
                axis equal
                axis tight
                %set(gca,'color',0*[1 1 1]);
                set(gca,'XTick',[],'YTick',[]);
                %ax.CLim = [0  maxVal];
                %colorbar('gray');

                title(t,sprintf('%s-MovBall-Positions.png',title_stim))
                %     xlabel(t,'X')
                %     ylabel(t,'Y');

                set(gcf,'PaperPositionMode','auto')
                
                fig = gcf;

                set(fig, 'Color', 'w');


                print(sprintf('%s-MovBall-Positions',title_stim),'-dpng');

                clearvars h
                close all
            end

            allM1(u,:,:) = nM1;
            allM2(u,:,:) = nM2;
            allM3(u,:,:) = nM3;
            allM4(u,:,:) = nM4;



            sumNeurons(1,u) = mean(nM1, "all");

            sumNeurons(2,u) = mean(nM3, 'all');

            sumNeurons(3,u) = mean(nM2, 'all');

            sumNeurons(4,u) = mean(nM4, 'all');

        end
       


    end

         %%% Normalize per direction
         offsetNames = unique(offsets);

        dirsSPKb = cell(length(offsetNames),length(direcNames));

        dirsSPKr = cell(length(offsetNames),length(direcNames));

        dirsSPKzb = cell(length(offsetNames),length(direcNames));

        dirsSPKza = cell(length(offsetNames),length(direcNames));

        minusRB = cell(length(offsetNames),length(direcNames));

        for d = 1:length(direcNames)

            for o = 1:length(offsetNames)

            stims = stimOn(directions==direcNames(d) & offsets==offsetNames(o));
                    
            %Baseline
            
            [Mb] = BuildBurstMatrix(goodU,round(p.t),round(stims-stimInter),round(stimInter)); 
            
            spkRateB = sum(Mb,3)/(stimInter/1000);

            dirsSPKb{o,d} = spkRateB;

            [MbAll] = BuildBurstMatrix(goodU,round(p.t),round(stimOn-stimInter),round(stimInter));

            spkRateBall = sum(Mb,3)/(stimInter/1000);

            %In order to deal with INF values in std and avoid INF
            %z-scores. 
            epsilon = 0.1; 
            
            %Response

            [Mr] = BuildBurstMatrix(goodU,round(p.t),round(stims),round(stimDur));
            
            spkRateR = sum(Mr,3)/(stimDur/1000);

            dirsSPKr{o,d} = spkRateR;

            %zScore
            AllSpk = (mean(spkRateB)+mean(spkRateR))./2;

            dirsSPKzb{o,d} = (spkRateR - mean(spkRateBall))./(std(spkRateBall)+epsilon);

              
            dirsSPKza{o,d} = (spkRateR - AllSpk)./(std(AllSpk)+epsilon);
            
            %R-B
            minusRB{o,d} = spkRateR - mean(spkRateB);

            end
            
        end           
  

    if heatmaps ==1
        animalSum{in_ind} = {allM1,allM2,allM3,allM4};
    end
        
        allDepth{in_ind} = {verticalDepth};

        allXdist{in_ind} = {XDist};
        
        DirSpkR{in_ind}=dirsSPKr;

        DirSpkB{in_ind}=dirsSPKb;

        DirSpkZB{in_ind}=dirsSPKzb;

        DirSpkZA{in_ind}=dirsSPKza;

        minusSpK{in_ind}=minusRB;

    in_ind = in_ind+1;
    end
end


if heatmaps ==1
    save(sprintf('%s-moving-ball-spkR',exp),"animalSum")
end
cd(basic_path+string(filesep)+exp)
save(sprintf('%s-unit-depths',exp),"allDepth")
save(sprintf('%s-unit-x-dist',exp),"allXdist")
save(sprintf('%s-DirSpkResponse',exp),"DirSpkR")
save(sprintf('%s-DirSpkBaseline',exp),"DirSpkB")
save(sprintf('%s-DirZScoreBaseline',exp),"DirSpkZB")
save(sprintf('%s-DirZScoreAll',exp),"DirSpkZA")
save(sprintf('%s-Good-Units',exp),"GoodUnits")
save(sprintf('%s-Res-Base',exp),"minusSpK")



%% Load results
cd(basic_path+string(filesep)+exp)


spkB = load(sprintf('%s-DirSpkBaseline.mat', exp));
spkR =load(sprintf('%s-DirSpkResponse.mat', exp));
spkZA =load(sprintf('%s-DirZScoreAll.mat', exp));
spkZB =load(sprintf('%s-DirZScoreBaseline.mat', exp));
allDepth =load(sprintf('%s-unit-depths.mat', exp));
goodUnits = load(sprintf('%s-Good-units.mat', exp));
RmB = load(sprintf('%s-Res-Base.mat', exp));




%% Selection of highest response values per offset

dirs = ["North","West","South","East"];

direcN = 4;

ofssetN = 9;

InN = 2; 

trials =15;

signU = cell(InN,direcN);

InsMax = cell(1,InN);

InsMaxO = cell(1,InN);

InsMaxB = cell(1,InN);

InsMaxA = cell(1,InN);

InsMinusA = cell(1,InN);

InsMinusB = cell(1,InN);

DepthMat = cell(1,InN);

InsMaxBt = cell(1,InN);

InsMinust = cell(1,InN);

InsMaxAt = cell(1,InN);

InsDivRB = cell(1,InN);

%Select top Minus R-B / Z-score baseline / All response off all offsets, per direction.

for i = 1:InN %insertions

    DepthMat{i} = cell2mat(allDepth.allDepth{i});

    maxValuesDIR =zeros(length(goodUnits.GoodUnits{i}), direcN);

    maxValuesDIRo = zeros(length(goodUnits.GoodUnits{i}), direcN);

    maxValuesZbDIR = zeros(length(goodUnits.GoodUnits{i}), direcN);

    maxValuesZaDIR =zeros(length(goodUnits.GoodUnits{i}), direcN);

    RmBbT_DIR = zeros(trials,length(goodUnits.GoodUnits{i}), direcN); %Response - baseline
    maxBZSTDIR = zeros(trials,length(goodUnits.GoodUnits{i}), direcN);
    maxAZSTDIR = zeros(trials,length(goodUnits.GoodUnits{i}), direcN); 

    RmBbM_DIR = zeros(length(goodUnits.GoodUnits{i}), direcN);

    divRB_DIR = zeros(length(goodUnits.GoodUnits{i}), direcN);

    for u = 1:length(goodUnits.GoodUnits{i})

        for d = 1:direcN
            
            maxVal = -inf;

            for o = 1:ofssetN

                %Selection of max z-score based on baseline:
                currentValueZB = mean(spkZB.DirSpkZB{i}{o,d}(:,u));
                
                currentValueM = mean(spkR.DirSpkR{i}{o,d}(:,u)) - mean(spkB.DirSpkB{i}{o,d}(:,u));

                %selection of max variable (mean response of trials per
                %offset)
                v = currentValueM;


                if v == inf %Sort out inf values.
                    v = 0;
                end

                if v > maxVal %select max value.
                    %baseline z-score mean max OR minus mean
                    maxVal= v;
                    %Maxval op offset
                    maxValO = mean(spkR.DirSpkR{i}{ofssetN-o+1,d}(:,u)) - mean(spkB.DirSpkB{i}{ofssetN-o+1,d}(:,u));

                    %baseline z-score mean
                    maxValzb = mean(spkZB.DirSpkZB{i}{o,d}(:,u));
                    %baseline Z-score per trial of max offset 
                    maxBZST = spkZB.DirSpkZB{i}{o,d}(:,u);
                    %minus trials
                    RmBbT = spkR.DirSpkR{i}{o,d}(:,u) - spkB.DirSpkB{i}{o,d}(:,u); %trials
                    %minus mean
                    RmBbM = mean(spkR.DirSpkR{i}{o,d}(:,u)) - mean(spkB.DirSpkB{i}{o,d}(:,u)); %mean
                    %div mean
                    DivRB = mean(spkR.DirSpkR{i}{o,d}(:,u))/mean(spkB.DirSpkB{i}{o,d}(:,u));            

                end

%                 %Selection of max z-score based on overall firing:
%                 currentValueZA = mean(spkZA.DirSpkZA{i}{o,d}(:,u));
% 
%                 if currentValueZA == inf %Sort out inf values.
% 
%                     currentValueZA = 0;
%                 end

%                 if currentValueZA > maxValza %select max value.
%                     %All z-score mean max
%                     maxValza = currentValueZA;
%                     %all Z-score per trial of max offset 
%                     maxAZST = spkZB.DirSpkZB{i}{o,d}(:,u);
% 
%                 end
%                 

                maxValuesDIR(u,d) = maxVal; %Max baseline Z-score across offset, per direction. 
                maxValuesDIRo(u,d) = maxValO;
                
                maxValuesZbDIR(u,d) = maxValzb;

                maxBZSTDIR(:,u,d) = maxBZST; 
                
                RmBbT_DIR(:,u,d) = RmBbT; % difference of spike rate for highest z-scores in offsets per trial
                RmBbM_DIR(u,d) = RmBbM; % difference of spike rate for highest z-scores in offsets

%                 maxValuesZaDIR(u,d) = maxValza; %Max overall response Z-score across offset, per direction. 
%                 maxAZSTDIR(:,u,d) = maxAZST;

                divRB_DIR(u,d) = DivRB;
                

            end
        end 
    end

    InsMax{i} = maxValuesDIR;

    InsMaxO{i} = maxValuesDIRo;

    InsMaxB{i} = maxValuesDIR;

    InsMinusB{i} = RmBbM_DIR;

    InsMaxBt{i} = maxBZSTDIR;

    InsMinust{i} = RmBbT_DIR;

    InsMaxA{i} =  maxValuesZaDIR;
    
    InsMaxAt{i} = maxAZSTDIR;

    InsDivRB{i} = divRB_DIR;
end

%% Plot figures with top Z-score per offset against depth
f = figure();

t = tiledlayout(InN,direcN,'TileSpacing','compact');

ts =1; 

maxZ = 0;

minZ = inf;

for i = 1:InN

    respU = sum(InsMaxB{i}>2.6,'all');

    for d = 1:direcN
        
        Val =  InsMax{i}(:,d);   

        ZB = InsMaxB{i}(:,d);

        ZA = InsMaxA{i}(:,d);

        MinusB = InsMinusB{i}(:,d);

        Div = InsDivRB{i}(:,d);

        yAxis = round(DepthMat{i})*-1;

        transp = 0.1;


        currentZM = max(ZB(~isinf(ZB))); 

        currentZm = min(ZB);

        if currentZM > maxZ
            maxZ = currentZM;
        end

        if currentZm < minZ
            minZ = currentZm;
        end

        tiles(ts) = nexttile;
        signif = ZB> 2.6; %99 % confidence


        scatter(MinusB(signif),yAxis(signif),[],ZA(signif), "filled", MarkerEdgeColor='k')
        hold on
        scatter(MinusB(~signif),yAxis(~signif), [], "filled", MarkerFaceAlpha=transp)
        grid on;
        set(gca, 'XGrid', 'off');
        ylabel(sprintf('In. %d: R.u/T.u = %d/%d',i,respU,length(goodUnits.GoodUnits{i})));
        yticks(-4000:500:0);
        xlabel(sprintf('R.u/T.u = %d/%d', sum(signif), length(goodUnits.GoodUnits{i})))
        hold off

        if ts ~= 1 & mod(ts-1, direcN)
              yticklabels("");
              yticks(-4000:500:0);
              set(gca,"YLabel",[]);
        end

        if ts <= direcN
            title(sprintf("%s",dirs(d)));
        end

        if ts <= direcN*InN - direcN
            xticks([])
        end

        signU{i,d} = goodUnits.GoodUnits{i}(signif);
        ts = ts+1;

    end
end

hcb = colorbar('Position', [0.94 0.11 0.02 0.815]); 
hcb.Title.String = "Z-score";% Adjust the position as needed
%hcb.Ticks([5:5:floor(maxZ/10)*10])

% Define the number of colors in the Sky colormap
numColors = 256;

% Create the custom Sky colormap
customSkyMap = zeros(numColors, 3);

% Define RGB values for the Sky colormap
lightBlue = [0.53, 0.81, 0.98]; % Light blue color
darkBlue = [0.0, 0.20, 0.40];   % Dark blue color

% Interpolate between light blue and dark blue
for i = 1:numColors
    customSkyMap(i, :) = (1 - (i - 1) / (numColors - 1)) * lightBlue + ((i - 1) / (numColors - 1)) * darkBlue;
end

colormap(customSkyMap);
caxis([2.6, ceil(maxZ)]);

title(t, sprintf('%s: Moving stimulus',animal));
xlabel(t,'Spike rate (spikes/sec) (response - baseline)');
ylabel(t,'Depth (um)');
linkaxes(tiles);
%linkprop(findall(gcf, 'Type', 'ColorBar'),'limits');
%linkprop(findall(gcf, 'Type', 'ColorBar'),'color');
prettify_plot;
screenSize = get(0, 'ScreenSize');
set(f, 'Position', screenSize);
print(f,sprintf('%s-Moving-ball',exp), '-dpng')
print(f,sprintf('%s-Moving-ball',exp), '-dpdf','-vector', '-fillpage')


%% Create colormaps

%figure;scatter(1:256,ones(1,256),[],(flipud(colormapR)+colormapB)/2,'filled')
%Create colormaps
white = [1.0, 1.0, 1.0];
black = [0, 0, 0];
dark_yellow = [0.9, 0.9, 0.0];
dark_red = [1, 0.0, 0.0];
dark_blue = [0.0, 0.0, 1];
dark_green = [0.0, 1, 0.0];

n_colors = 64;  % Number of colors in the colormap
colormapG = zeros(n_colors, 3);
colormapB = zeros(n_colors, 3);
colormapY = zeros(n_colors, 3);
colormapR = zeros(n_colors, 3);


for i = 1:n_colors
    colormapG(i, :) = (1 - (i - 1) / (n_colors - 1)) * black + ((i - 1) / (n_colors - 1)) * dark_green;
    colormapB(i, :) = (1 - (i - 1) / (n_colors - 1)) * black + ((i - 1) / (n_colors - 1)) * dark_blue;
    colormapY(i, :) = (1 - (i - 1) / (n_colors - 1)) * black + ((i - 1) / (n_colors - 1)) * dark_yellow;
    colormapR(i, :) = (1 - (i - 1) / (n_colors - 1)) * black + ((i - 1) / (n_colors - 1)) * dark_red;
end

RB = (flipud(colormapR)+colormapB);

BG = (flipud(colormapB)+colormapG);

GY = (flipud(colormapG)+colormapY);

YR = (flipud(colormapY)+colormapR);



color_triplets = [RB; BG; GY; YR]; % Replace with your desired color triplet

color_tripletsW = zeros(length(color_triplets),n_colors*3);

[allC wRange] = size(color_tripletsW);


%% Calculate colormaps

[CM,CMP,h, cmapR_theta]=colormap2D('plotPolarExample',1);

%prettify_plot

figure(7)



%% Calculate strength of selectivity. (Color) (N-S) = y, (W-E) = x

%R = squeeze(mean(InsMinust{1}));

NeuronColor = cell(1,InN);

[nRC,~,nTC]=size(cmapR_theta);

th = linspace(0,2*pi,50);
r = 2; %signif;

%t = tiledlayout(1,InN);

for in = 1:InN
    

    NeuronColor{in} = zeros(length(goodUnits.GoodUnits{in}),3);
    R = InsMax{in}; 
    
    Ro = InsMaxO{in};

    %Find closest points to polar coordinates given the colorwheel script and

    %select = [(R(:,2)-Ro(:,4))-abs(Ro(:,2)-R(:,4))  (R(:,1)-Ro(:,3))-abs(Ro(:,1)-R(:,3))];

    select = [(R(:,2)-R(:,4))  R(:,1)-R(:,3)];

    %selectN = [xNorm,yNorm];

    [thetaR, rhoR] = cart2pol(select(:,1), select(:,2));

    
    divFactor = 2*pi;

    thetaR(thetaR<0) = abs(thetaR(thetaR<0))+pi;

    maxR =8;
    rhoR(rhoR>maxR) = maxR;

    hf = figure('Position',[400*in 400 400 450]);


    cMapRed = zeros(nR,3,20);
    NumU =length(thetaR);

    rhoR(rhoR == 0) = 0.001; 
    
    thetaR(thetaR == 0) = 0.001;%prevents index from being 0

    hP=polaraxes;colormap(hP,'hsv');hold on;

    
    selU = sum(rhoR>2);
  

    %Selective units
    ns = find(round(rhoR*10)/10==8.3);
    rad2deg(thetaR(ns))
    rhoR(59)

    rho(59)

    for i= 1:NumU

            
        polarplot(thetaR(i),rhoR(i),'.','MarkerSize',20, ...
            'Color', cmapR_theta( ceil((thetaR(i)/(divFactor))*nRC), : , ceil((rhoR(i))/maxR*nTC)) )

%         p = polarplot(th,r+zeros(size(th)),'-','LineWidth',3);
%         p.Color =[0.7  0.7  0.7 0.2];

        NeuronColor{in}(i,:) = cmapR_theta( ceil((thetaR(i)/(divFactor))*nRC), : , ceil((rhoR(i))/maxR*nTC) ) ;
        
    end
    % Get the current figure and axes
    fig = gcf;
    ax = gca;

    % Set the color of the figure and axes to black
    set(fig, 'Color', 'k');
    set(ax, 'Color', 'K');
   

% Change the color of the text and grid to white
    ax.RColor = [1 1 1];
    ax.ThetaColor = [1 1 1];
    ax.MinorGridColor= [1 1 1];
    ax.GridColor = [1 1 1];
    ax.MinorGridAlpha = 0.8;
    ax.GridAlpha = 0.8;
    ax.FontSize = 14;

    %r lims
    rlim([0 maxR]);
    title(sprintf('In. %d: S.u/T.u = %d/%d',in, selU,NumU),'Color','w')
    hold off
   % prettify_plot;

end





%%





%Find closes value

for u = 1:length(tethaR)

    [minT, Nind(u)] = min(abs(thetaO-tethaR(u)));

end


%X

rRN = normalize(rhoR,'range');

figure
%colormap(cmapa);
for u = 1:length(Nind)
    
    ps = polarscatter(thetaR(u),rhoR(u),'filled');
    ps.MarkerFaceColor = cmapaN(Nind(u),:);
    ps.MarkerFaceAlpha =rRN(u);
    hold on
end
hold off


figure
scatter(invUx(unitIx),invUy(unitIy(u)),[],neuronColor','filled')

figure
imagesc(rgb)
axis equal

color =squeeze(rgb(401,401,:))';
figure
scatter(401,401,[],color', 'filled')

%%

squeeze(rgb(:,:,1));



%%

%%

for c = 1:allC


 color_tripletsW(c,:) = reshape([linspace(color_triplets(c,1), 1, n_colors)', ...
                       linspace(color_triplets(c,2), 1, n_colors)', ...
                       linspace(color_triplets(c,3), 1, n_colors)']',1,[]);  
end

colorCircle = cell(128,128);

%RB

for i = 1:n_colors

    for j = 1:n_colors
        for c = 1:n_colors
            for w = 1:3:wRange
                colorCircle{n_colors+i,n_colors}  =  color_tripletsW(c,w:w+2);
            end
        end
    end

end


%cartesian matrix.

x = repmat(linspace(-1,1,n_colors*2),n_colors*2,1);
y = repmat(linspace(-1,1,n_colors*2)',1,n_colors*2);


% Normalize the Cartesian coordinates to have a minimum of -1 and a maximum of 1
x_normalized = x / ((matrix_size-1)/2);
y_normalized = y / ((matrix_size-1)/2);

% Create the Cartesian coordinate matrix
cartesian_matrix = sqrt(x_normalized.^2 + y_normalized.^2);

% Display the result
imshow(cartesian_matrix, []);
title('Normalized Cartesian Coordinate System (-1 to 1)');colormap(colormapR+colormapB);
colorbar
%







%%

unitsI = zeros(1,InN);
for in = 1:InN
    unitsI(in) = length(goodUnits.GoodUnits{in});
end

maxUnits = max(unitsI); 

figure()
T = tiledlayout(direcN,1,'TileSpacing', 'compact');

for d = 1:direcN

[m, ind] = max(squeeze(mean(InsMaxBt{1},1)),[],2);

color_names = cell(size(ind));

for i = 1:numel(ind)
    switch ind(i)
        case 1 %North
            color_names{i} = colormapB; %[0 0.4470 0.7410];
        case 2 %West
            color_names{i} = colormapR; %[0.4660 0.6740 0.1880];
        case 3 %South
            color_names{i} = colormapG; %[0.6350 0.0780 0.1840];
        case 4 %East
            color_names{i} = colormapY; %[0.9290 0.6940 0.1250];
        otherwise
            color_names{i} = colormapY; % Optional: Handle other values if needed
    end
end

[y x d] = size(InsMaxBt{in});

[si indS] = sort(ind,'descend');


meanZf = InsMaxBt{in};

limitC = 10;
meanZf(meanZf>limitC) = limitC; 
meanZf(meanZf<2) = 0; 

meanZFbig = zeros(y,maxUnits,4);

meanZFbig(:,1:x,:) = meanZf;

meanZf = meanZFbig;

SmeanZf = meanZf(:,indS,:);

[y x d] = size(meanZf);

nexttile(T) 
t = tiledlayout(1,x, 'Parent', T, 'TileSpacing', 'none');
t.Layout.Tile = in;
%t.Layout.TileSpan = [1 1];
checkTiles = zeros(1,4);
tile = cell(1,numel(indS));

si(si == 0) = 1;

for c = 1:numel(indS)
    %xline(c,'Color',color_names{c},'Alpha', 0.3, 'LineWidth',5)
    tile{c} = nexttile(t);
    
    imagesc(1,1:y,meanZf(:,indS(c),si(c)));
    colormap(tile{c},color_names{indS(c)});
    caxis([0, limitC]);
    set(gca,'Visible','off');

%     if c ~= 1
%         yticks(tile{c},[]);
%         
%     end
%     if c<numel(indS) 
%         xline(numel(indS)+0.5)
%     end

    xticks([]);    
    
    
end
%title('Blue = North; West = green; South = Red; East = Yellow')
ylabel(t, sprintf('Insertion %d',in))
xticks([1:x]);



% caxis([0, 10]);
% c.Title.String = ('z-score');
end

prettify_plot;
ylabel(t,'15 Trials')
xlabel(t, 'Neurons ordered by direction selectivity')




%%

spkR.DirSpkR{1}{1}(goodUnits.GoodUnits{1}==292)

spkZB.DirSpkZB{1}{1}(goodUnits.GoodUnits{1}==292)

goodUnits.GoodUnits{1}(spkZB.DirSpkZB{1}{1}>2 & spkZB.DirSpkZB{1}{1}<inf)


%%

%List of responsive neurons
%PV67_1 = [49, 70, 99, 132, 136, 151, 158, 166, 174, 177, 205, 239, 240, 246, 254, 265, 267, 272, 287, 298, 301, 305, 306, 336, 345, ]

%Summary accordying to depth
animalSum = load('PV102_experiment_18_7_23-moving-ball-spkR.mat');

animalSum = animalSum.animalSum;

numDir = 4;

figure()
tiledlayout(length(ins), numDir,'TileSpacing','tight')

allDepth = load('PV102_experiment_18_7_23-unit-depths.mat');

direcs = ["North","West","South","East"];

for i = ins %insertions

    
    insertion = mean(cell2mat(animalSum{i}),[2 3],"omitnan");
    AllN =[];
    colors = rand(length(allDepth{i}{1}),3);
    %colors(colors(:, 1) < 0.2, :) = [1, 0, 0];
    for d = 1:numDir %directions
        
        zscoreMat = animalSum{i}{d};

        meanSCu = zeros(1,length(zscoreMat)); 
        
        for u = 1:length(zscoreMat)
            uZS = zscoreMat(u,:,:);
            
            NuZS = (uZS-mean(uZS,"all"))./std(uZS); %Double normalization.

            if sum(uZS > 10,'all') > 0
                meanSCu(u) = mean(uZS(uZS > 10),'all');
            else
                meanSCu(u) = mean(uZS, "all");
            end
        end

        %AllN = [AllN;mean(zscoreMat(zscoreMat>2),[2 3],"omitnan")'];

        AllN = [AllN;meanSCu];


        AllNnorm = normalize(AllN,"range");

        %xAxis = rand(1,length(verticalDepth));
        yAxis = round(allDepth{i}{1})*-1;
       

        szFactor = 200;
        transp = 0.1;
        hold on
        

        nexttile
        zscorePerU = AllN(d,:);
        signif = zscorePerU > 2;
        
        scatter(zscorePerU(signif),yAxis(signif), [], colors(signif), "filled")
        hold on 
        title(sprintf("%s: Insertion %d - direction %s",animal,i, direcs(d)))
        scatter(zscorePerU(~signif),yAxis(~signif), [], colors(~signif), "filled",MarkerFaceAlpha=transp)
        %xlim([min(insertion,[], 'all')-1 max(insertion,[], 'all')+1])
        hold off
% scatter(xAxis,yAxis,round(abs(AllN(1,:))*szFactor)+1,"filled","MarkerFaceColor",[0 0.4470 0.7410])
% scatter(xAxis+1.5,yAxis,round(abs(AllN(3,:))*szFactor)+1,"filled","MarkerFaceColor",[0.3010 0.7450 0.9330])
% scatter(xAxis+3,yAxis,round(abs(AllN(2,:))*szFactor)+1,"filled","MarkerFaceColor",[0.9290 0.6940 0.1250])
% scatter(xAxis+4.5,yAxis,round(abs(AllN(4,:))*szFactor)+1,"filled","MarkerFaceColor",[0.8500 0.3250 0.0980])
    end
    
end
% nexttile
% scatter(AllN(3,:),yAxis,"filled", "MarkerFaceColor",[0.4660 0.6740 0.1880])
% nexttile
% scatter(AllN(2,:),yAxis,"filled", "MarkerFaceColor",[0.6350 0.0780 0.1840])
% nexttile
% scatter(AllN(4,:),yAxis,"filled", "MarkerFaceColor",[0.9290 0.6940 0.1250])




% scatter(x_pos,yAxis,round(abs(AllN(1,:))*szFactor)+1,"filled","MarkerFaceColor",[0 0.4470 0.7410],MarkerFaceAlpha=transp)
% scatter(xAxis,yAxis,round(abs(AllN(3,:))*szFactor)+1,"filled","MarkerFaceColor",[0.4660 0.6740 0.1880],MarkerFaceAlpha=transp)
% scatter(xAxis+1.1,yAxis,round(abs(AllN(2,:))*szFactor)+1,"filled","MarkerFaceColor",[0.6350 0.0780 0.1840],MarkerFaceAlpha=transp)
% scatter(xAxis+1.1,yAxis,round(abs(AllN(4,:))*szFactor)+1,"filled","MarkerFaceColor",[0.9290 0.6940 0.1250],MarkerFaceAlpha=transp)
%xlim([-0.5 6])
ylim([-ceil(InserDepth/1000)*1000,0])
%xticks([0.5 2 3.5 5])
ylabel("depth (um)")
%xticklabels(["North","South","West","East"]);
legend(["North","South","West","East"],'location','southoutside')
yline(-(InserDepth-3840),'-b', 'Start of active Chs', 'LabelHorizontalAlignment', 'center', LabelVerticalAlignment='middle', LineWidth=1.5, HandleVisibility='off')
set(gca, 'YGrid', 'on', 'XGrid', 'off')

print(sprintf("Resp-depth-%s-Insertion %d",animal,in),'-dpdf','-vector');
print(sprintf("Resp-depth-%s-Insertion %d",animal,in),'-dpng');
%yticklabels(sort(0:400:ceil(InserDepth/1000)*1000,"descend"))

%%


%Summary off all neurons:

Northmean = (squeeze(mean(allM1,1,"omitnan")));
Westmean = (squeeze(mean(allM2,1,"omitnan")));
Southmean = (squeeze(mean(allM3,1,"omitnan")));
Eastmean = (squeeze(mean(allM4,1,"omitnan")));



concatMean = cat(2, Northmean,Westmean,Southmean,Eastmean);


t = tiledlayout(2, 2,'TileSpacing','tight');

LB=flipud(lbmap(256,'BrownBlue'));

axisL = 1:length(nM1);

nexttile
imagesc(Northmean);
%plot it
title('North');
colormap(LB);
caxis([-max(abs(concatMean(:))) max(abs(concatMean(:)))]);
ax = gca;
set(gca,'XTick',[],'YTick',[]);
axis equal
axis tight

nexttile

imagesc(Southmean);
%plot it
title('South');
colormap(LB);
caxis([-max(abs(concatMean(:))) max(abs(concatMean(:)))]);
ax = gca;
set(gca,'XTick',[],'YTick',[]);
% ax.CLim = [0  maxVal];
axis equal
axis tight
colorbar();


nexttile
imagesc(Eastmean); %plot it
title('East');
colormap(LB);
caxis([-max(abs(concatMean(:))) max(abs(concatMean(:)))]);
ax = gca;
set(gca,'XTick',[],'YTick',[]);
axis equal
axis tight
%ax.CLim = [0  maxVal];
%colorbar('gray');


nexttile
imagesc(Westmean); %plot it
title('West');
colormap(LB);
caxis([-max(abs(concatMean(:))) max(abs(concatMean(:)))]);
ax = gca;
set(gca,'XTick',[],'YTick',[]);
axis equal
axis tight
%ax.CLim = [0  maxVal];
%colorbar('gray');

title(t,sprintf('%s-MovBall-Positions-mean',title_stim))
%     xlabel(t,'X')
%     ylabel(t,'Y');

%set(gcf,'PaperPositionMode','auto')
print(sprintf('%s-MovBall-Positions-mean',title_stim),'-dpdf','-vector');


%%


matAll1 = mean(allM(1,:));
CN = normalize(sumNeurons);

N = normalize(CN,2);

% set(gca,'xticklabel',names)
%yticklabels(['North';'South';'West';'East']);
figure
imagesc(Nrthmean);
title(sprintf('%s-MovBall-Positions-all-neurons',title_stim))
colorbar
set(gca, 'Ytick' ,[1;2;3;4], 'YTickLabel',{'North';'South';'West';'East'});
print(sprintf('%s-MovBall-Positions-all-neurons.png',title_stim),'-dpng');

%
%%
function r=toWhite(color_triplet)

    r = [linspace(color_triplet(1), 1, num_colors)', ...
                       linspace(color_triplet(2), 1, num_colors)', ...
                       linspace(color_triplet(3), 1, num_colors)'];  
end
