

% path = '\\sil3\data\\Large_scale_mapping_NP\lizards\PV139\PV139_Experiment_6_2_24\Insertion1\catgt_PV139_Experiment_6_2_24_1_g0';
%
% NP = NPAPRecording(path);

%%%%%Function inputs = NP, stimName, ismoving;

%stimName = "MovBall";

%stimName = "rectGrid";

function [onsetSync offsetSync onSync offSync] = NPdiodeExtract(NP,newDiode,ismoving,stimName,ttl_index,digCH,syncCH) %dig -16/1 sync 0/6

%
%%%%MATLAB PTB SPECS: %%%%%%%%%%%%%%%%%%%%%%
cd(NP.recordingDir)

validOptions = {'MB', 'SDG', 'BB', 'FFF','RGN','RG','NS','NSC'};

% Validate the input
stimName = validatestring(stimName, validOptions);

% Example logic based on the validated input
switch stimName
    case 'MB' %%Linearly moving ball
        stimtype= 'linearlyMovingBall';

    case 'SDG' %%Static and drifting gratings
        stimtype = 'StaticDrifting';

    case 'BB' %%Linearly moving bouncing balls
        stimtype = 'linearlyMovingBouncing';

    case 'RG' %%Rectangle grid
        stimtype = 'rectGrid';

    case 'FFF' %%Full field flash
        stimtype = 'fullFieldFlash';

    case 'RGN' %%Rectangle Grid Noise
        stimtype = 'rectNoiseGrid';

    case 'NS' %%Novelty stim
        stimtype = 'NoveltyStim';

    case 'NSC' %%Novelty stim control
        stimtype = 'NoveltyStimC';

    otherwise
        error('Invalid input option.');
end

fileName = "DIODE_"+string(stimName);

MovExist =  exist(NP.recordingDir+"\"+fileName+"OnFrame"+".txt",'file') & exist(NP.recordingDir+"\"+fileName+"OffFrame"+".txt",'file')...
        & exist(NP.recordingDir+"\"+fileName+"Onset"+".txt",'file') & exist(NP.recordingDir+"\"+fileName+"Offset"+".txt",'file') & ismoving;

nonMovExist = exist(NP.recordingDir+"\"+fileName+"Onset"+".txt",'file') & exist(NP.recordingDir+"\"+fileName+"Offset"+".txt",'file');

if MovExist && newDiode == 0

    onsetSync = readmatrix(fileName+"Onset"+".txt");
    offsetSync =  readmatrix(fileName+"Offset"+".txt");
    onSync = readmatrix(fileName+"OnFrame"+".txt");
    offSync = readmatrix(fileName+"OffFrame"+".txt");

elseif nonMovExist && newDiode == 0

    onsetSync = readmatrix(fileName+"Onset"+".txt");
    offsetSync =  readmatrix(fileName+"Offset"+".txt");

else


    patternIndex = strfind(string(NP.recordingDir), "\catgt");

    endIndex = patternIndex(1)-1;
    stimDir = string(NP.recordingDir);
    stimDir = extractBetween(stimDir,1,endIndex);

    file = dir (stimDir);
    filenames = {file.name};


    %%% add SDG and BB

    % Define the valid options
    

    stimFiles = filenames(contains(filenames,stimtype));
    j =1;
    if size(stimFiles) ~= [0 0]
        for i = stimFiles
            stim = load(stimDir+"\"+string(i));
            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    %%%%%GET UN-SYNCED DIGITAL TRIGGERS %%%%%%%%%%%%%%%%%%%%%%


    file = dir (NP.recordingDir);
    filenames = {file.name};

    cd(NP.recordingDir)

    Sstart = readmatrix(string(filenames(contains(filenames,sprintf("_tcat.nidq.xd_%d_1_0",digCH)))))*1000;
    Send = readmatrix(string(filenames(contains(filenames,sprintf("_tcat.nidq.xid_%d_1_0",digCH)))))*1000;
    if Send(1)<Sstart(1)
        Send = Send(2:end);
    end
    onset = readmatrix(string(filenames(contains(filenames,sprintf("_tcat.nidq.xd_%d_2_0",digCH)))))*1000;
    offset = readmatrix(string(filenames(contains(filenames,sprintf("_tcat.nidq.xid_%d_2_0",digCH)))))*1000;




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


    stimInter = mean(stimOn(2:end)-stimOff(1:end-1));

    stimOn = stimOn';

    stimDur = mean(stimOff- stimOn');%transform into horizontal vector.


    %%%%%MARK'S FUNCTION (NOT USED)

    % [frameShifts,upCross,downCross,diffStats,transitionNotFound,T] = frameTimeFromDiode(NP,'tStart',stimOn(1),'tEnd',stimOff(end),'analogChNum',1,'noisyAnalog',true);
    %
    % %%Stim onset and offset from analog signal
    % MBcrossUp = upCross(upCross >stimOn(1)-400 & upCross<stimOff(end)+400);
    % MBcrossDown = downCross(downCross >stimOff(1)-400 & downCross<stimOff(end)+400);
    %%%%%%%

    %%%%%%GET ANALOG DATA:%%%%%%%%%%%%%%%%%%%%%%
    startS = stimOn(1)-50;

    All = NP.getAnalogData(1,startS,stimOff(end)-startS+100);

    AllD = squeeze(All);

    frameRate = 60;

    frameSamples=round(NP.samplingFrequencyNI/frameRate);


    %%%% IF MOVING CALCULATE FRAMES %%%%%%%%%%%%%%%%%%%%%%

    if ismoving == 1
        onFrame = [];
        offFrame =[];
        onset = [];
        offset = [];

%         if string(stimName) == "MB"
%             nFrames = unique(cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'nFrames')))));
        if string(stimName) == "SDG"
             nFrames = size(cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'flipOnsetTimeStamp')))),2);
        end
        
        FramesPsec = cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'fps'))));%ball.VSMetaData.allPropVal{60,1};

        samplesPframe = 1/FramesPsec*(NP.samplingFrequencyNI);


        if string(stimName) == "SDG"

            staticTime = cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'static_time'))));
        end


        for i = 1:length(stimOn)

            %%%%Initialize variables and define start and end of snip of
            %%%%signal:

            pklocUpN = [];
            pklocDownN = [];
            %i=44;
            startSnip  = round((stimOn(i)-stimOn(1))*(NP.samplingFrequencyNI/1000))+1;
            endSnip  = round((stimOff(i)-stimOn(1)+1000)*(NP.samplingFrequencyNI/1000));

            interstim = cell2mat(stim.VSMetaData.allPropVal(find(strcmp(stim.VSMetaData.allPropName,'interTrialDelay'))));

            if interstim <1 %Check if interstim is less than 1 and modify endSnip
                endSnip  = round((stimOff(i)-stimOn(1)+200)*(NP.samplingFrequencyNI/1000));
            end

             %%%%Apply median filter and make sure that selection does not go out of bounds:
            if endSnip>length(AllD)
                signal = AllD(startSnip:end);
                fDat=medfilt1(signal,15);      
            else
                signal =(AllD(startSnip:endSnip));
                fDat=medfilt1(signal(1:end-1),15); %median dilter is used to mantained edges and facilitate downshifts more efectively         
            end

            %%%%Apply low pass filter:

            % Use the designfilt function to create a low-pass filter
            d = designfilt('lowpassiir', 'FilterOrder', 4, ...
                'HalfPowerFrequency', frameRate+5, 'SampleRate', NP.samplingFrequencyNI);
            % Apply the filter using filtfilt to preserve phase
            
            fDat = filtfilt(d, fDat);

            %%%%Add higher or lower values at the end of the signal for last frame to be detected as
            %%peak:

            n=10; %samples to add at the end for last signal be recognized as peaks
            y = 0.01;%Step in which additional samples change value

            if mean(signal(end)-1000)> signal(end)

                signal = [signal;signal(end) + (1:n-1)' * y];
                fDat = [fDat;fDat(end) + (1:n-1)' * y];
            
            else %mean(signal(end)-1000)< signal(end)

                signal = [signal;signal(end) -  (1:n-1)' * y];
                fDat = [fDat;fDat(end) - (1:n-1)' * y];
            end

            %%%%Detect peaks
            stdS = std(signal);
            [pkvalsUp, pklocUp] = findpeaks(fDat,'MinPeakProminence',0.8*stdS,'MinPeakDistance',samplesPframe*2-samplesPframe*0.3);

            [pkvalsDown, pklocDown] = findpeaks(fDat*-1,'MinPeakProminence',0.8*stdS,'MinPeakDistance',samplesPframe*2-samplesPframe*0.3);

            %%%%filter out outlier peaks If there are more peaks than max
%             
%             tic
            numSegments = 2;


            % Initialize array to store change points from each segment
            change_points = cell(1, numSegments);
            %
            %             % Use parfor for parallel processing
            dFactor =10;
            dsFDat = downsample(fDat,dFactor);
            segmentLength = length(dsFDat) / numSegments;
            %
            tic
            j = 0;
            for par = 1:numSegments
                segmentData = dsFDat((par-1)*segmentLength + 1:par*segmentLength);
                [change_point ~] = findchangepts(single(segmentData), 'MaxNumChanges', 1, 'Statistic', 'std');
                if ~isempty(change_point)
                    change_points{par} = (change_point+length(segmentData)*j*par)*dFactor;
                end
                j=j+1;
            end
            toc

            %             % Combine results from all segments
            %             try
            %                 cpts = cell2mat(change_points) + repelem(0:segmentLength:(numSegments-1)*segmentLength, cellfun(@numel, change_points));
            %             catch
            if isempty(change_points{2})
                change_points{2} =  fDat(end);
            end
            if isempty(change_points{1})
                change_points{1} = 1;
            end
            cpts = cell2mat(change_points);
%             end
            %toc
            %[cpts1, ~] = findchangepts(single(fDat), 'MaxNumChanges', 2, 'Statistic', 'std'); %Detect signifficant changes in signal statistics
            %3 main changes happen in SDG, the On signal, the start of frames, and
            %the off signal. 
            toc
            if string(stimName) == "MB"
                pklocUp = pklocUp(pklocUp<max(cpts));

                pklocDown = pklocDown(pklocDown<max(cpts));

                pklocUp = pklocUp(pklocUp>min(cpts));

                pklocDown = pklocDown(pklocDown>min(cpts));
            end

            %%%%To eliminate risk of noise peaks in SDG, allow frame shifts
            %%%%only after static phase of gratings
            if string(stimName) == "SDG"
                 [cpts, ~] = findchangepts(fDat, 'MaxNumChanges', 4, 'Statistic', 'mean'); 

                try
                    if length(pklocUp) > nFrames/2
                        pklocUpM = pklocUp(pklocUp>cpts(2)&pklocUp<=cpts(3));
                        pklocDownM = pklocDown(pklocDown>cpts(2)&pklocDown<=cpts(3));

                        if length(pklocUpM)<nFrames/3
                            pklocUpM = pklocUp(pklocUp>cpts(1)&pklocUp<=cpts(2));
                        end
                        if length(pklocUpM)<nFrames/3
                            pklocUpM = pklocUp(pklocUp>cpts(2)&pklocUp<=cpts(4));
                        end
                        if length(pklocDownM)<nFrames/3
                            pklocDownM = pklocDown(pklocDown>cpts(1)&pklocDown<=cpts(2));
                        end
                        if length(pklocDownM)<nFrames/3
                            pklocDownM = pklocDown(pklocDown>cpts(2)&pklocDown<=cpts(4));
                        end


                        if mean(signal(1:cpts(1))) > mean(signal(cpts(1):cpts(1)+100-n)) %signal starts in down shift
                            try
                                pklocDown = [min(cpts(1),pklocDownM(1)-staticTime*NP.samplingFrequencyNI); pklocDownM];
                            catch
                                2+2
                            end
                            pklocUp = pklocUpM;
                        else
                            %pklocDown = [pklocDownM;cpts(3)];%signal starts in up shift
                            pklocUp = [min(cpts(1),pklocUpM(1)-staticTime*NP.samplingFrequencyNI);pklocUpM];
                            pklocDown = pklocDownM;
                        end
                    end
                catch
                    2+2
                end
            end

            if isempty(pklocUpN)
                pklocUpN = pklocUp;
            end
            if isempty(pklocDownN)
                pklocDownN = pklocDown;
            end

            %%%% Assign trial On and Off . 

            try

            [startP indS] =  min([pklocDownN(1),pklocUpN(1)]);

            catch

                [pkvalsUp, pklocUp] = findpeaks(fDat,'MinPeakProminence',0.8*stdS,'MinPeakDistance',samplesPframe*2-samplesPframe*0.3);

                [pkvalsDown, pklocDown] = findpeaks(fDat*-1,'MinPeakProminence',0.8*stdS,'MinPeakDistance',samplesPframe*2-samplesPframe*0.3);

       %%%Quick and dirty solution for failed trials=few
                pklocUpN = pklocUp;
                pklocDownN = pklocDown;

                [startP indS] =  min([pklocDownN(1),pklocUpN(1)]);


            end

            if indS ==1
                startP = pklocDownN;
            else
                startP = pklocUpN;
            end

            [endP indE] =  max([pklocDownN(end),pklocUpN(end)]);

            if indE ==1
                endP = pklocDownN;
            else
                endP = pklocUpN;
            end


            offFrame = [offFrame,endP'./(NP.samplingFrequencyNI/1000)+stimOn(i)];

            offset = [offset endP(end)./(NP.samplingFrequencyNI/1000)+stimOn(i)];

            onFrame = [onFrame, startP'./(NP.samplingFrequencyNI/1000)+stimOn(i)];

            onset = [onset startP(1)./(NP.samplingFrequencyNI/1000)+stimOn(i)];

            if endP(end)./(NP.samplingFrequencyNI/1000)-startP(1)./(NP.samplingFrequencyNI/1000) <2000

                2+2

            end



        end
%%%%Test:

%%% Check other digital triggers:

% figure()
% plot(segmentData);xline(change_point);
% 
% %(endP(end)-startP(1))/NP.samplingFrequencyNI
% figure()
% plot(signal);hold on; plot(fDat)
% 
% %[cpts, ~] = findchangepts(signal, 'MaxNumChanges', 2, 'Statistic', 'std');
% % xline(startP(1))
% %     xline(endP(end))
% xline(cpts);
% xline([pklocDown'])
% xline([pklocUp'],'red')
% xline(change_points{1})%*dFactor)
%hold on; plot(round((TTL-stimOn(i))*(NP.samplingFrequencyNI/1000)),ones(size(TTL))*max(fDat),'ro', 'MarkerSize', 0.2, 'LineWidth', 2,'MarkerFaceColor','b')


    else %Notmoving

        onset = [];
        offset = [];
        missedStims = [];
        for i =1:length(stimOn)

            startSnip  = round((stimOn(i)-stimOn(1)+10)*(NP.samplingFrequencyNI/1000))+1;
            startSnip(startSnip<1) = 1;
            endSnip  = round((stimOff(i)-stimOn(1)+500)*(NP.samplingFrequencyNI/1000));

            if endSnip>length(AllD)
                signal = AllD(startSnip:end);
                fDat=medfilt1(signal,2*frameSamples);
                flipD = flip(fDat);
            else
                signal =(AllD(startSnip:endSnip));
                fDat=medfilt1(signal,2*frameSamples);
                flipD = flip(fDat);
            end


%             try
            %stdS = std(fDat);
            stdD = std(diff(signal));
            [pkvalOn, pklocOn] = findpeaks(diff(-1*fDat),'MinPeakProminence',stdD,'MinPeakDistance',(stimDur-300)*NP.samplingFrequencyNI/1000);

            [pkvalOff, pklocOff] = findpeaks(diff(-1*flipD),'MinPeakProminence',stdD,'MinPeakDistance',(stimDur-300)*NP.samplingFrequencyNI/1000);

            pklocOn = min(pklocOn);

            posInverted =[1:length(fDat);flip(1:length(fDat))];
            pklocOff = posInverted(2,min(pklocOff));

            if isempty(pkvalOn)
                missedStims = [missedStims i];
                pklocOn = stimOn(i)*(NP.samplingFrequencyNI/1000);
                pklocOff = stimOff(i)*(NP.samplingFrequencyNI/1000);

                onset = [onset pklocOn./(NP.samplingFrequencyNI/1000)];
                offset = [offset pklocOff./(NP.samplingFrequencyNI/1000)];
            else

                onset = [onset pklocOn./(NP.samplingFrequencyNI/1000)+stimOn(i)];
                offset = [offset pklocOff./(NP.samplingFrequencyNI/1000)+stimOn(i)];
            end
%             afterOn = find(signal > mean(pkvalOn) - 5*std(pkvalOn));
%             pklocOfflist = find(afterOn>pklocOn+100*NP.samplingFrequencyNI/1000);
%             pklocOff = afterOn(pklocOfflist(1));
              
%             catch
%                 print(i)
%             end

        end

%             figure()
%             plot(signal)
%             hold on;plot(fDat)
%             $hold on;plot(diff(fDat))
%             xline(pklocOn,'red')
%             xline(pklocOff)
%             yline(mean(fDat))
%             xline((stimOn(1:5)-stimOn(1))*NP.samplingFrequencyNI/1000)

    cd(NP.recordingDir)
    fileID = fopen(fileName+"missedStims"+".txt",'w');
    fprintf(fileID, '%d\n', missedStims);
    fclose(fileID);

    end



    %%%% SYNC %%%%%%%%%%%%%%%%%%%%%%

    if nargin>=7 %% If there is a sync waveform:

        % Search in for neural SQW and NI SQW and load them

        Neur = readmatrix(dir(fullfile(NP.recordingDir, '*imec0.ap.xd_384_6_500.txt*')).name);
        originSQW = readmatrix(dir(fullfile(NP.recordingDir, sprintf('*nidq.xd_%d_%d_500.txt*',digCH,syncCH))).name);
        if ismoving
            offSync = interp1(originSQW'*1000, Neur'*1000, offFrame, 'linear');
            onSync = interp1(originSQW'*1000, Neur'*1000, onFrame, 'linear');

            fileID = fopen(fileName+"OnFrame"+".txt",'w');
            fprintf(fileID, '%d\n', onSync);
            fclose(fileID);

            fileID = fopen(fileName+"OffFrame"+".txt",'w');
            fprintf(fileID, '%d\n', offSync);
            fclose(fileID);

        end
        onsetSync = interp1(originSQW'*1000, Neur'*1000, onset, 'linear');
        fileID = fopen(fileName+"Onset"+".txt",'w');
        fprintf(fileID, '%d\n', onsetSync);
        fclose(fileID);

        offsetSync = interp1(originSQW'*1000, Neur'*1000, offset, 'linear');
        fileID = fopen(fileName+"Offset"+".txt",'w');
        fprintf(fileID, '%d\n', offsetSync);
        fclose(fileID);



        % %test sync
        %
        % test = readmatrix('PV139_Experiment_6_2_24_1_g0_tcat.nidq.xd_1_2_0.txt');
        % testS = interp1(originSQW, Neur, test, 'linear');
        %
        % testTprime = readmatrix('out_2.txt');
        %
        % t = [testTprime,testS];
        %
        % t(end,:)*1000

    else

        fileID = fopen(fileName+"Onset"+".txt",'w');
        fprintf(fileID, '%d\n', onset);
        fclose(fileID);

        fileID = fopen(fileName+"Offset"+".txt",'w');
        fprintf(fileID, '%d\n', offset);
        fclose(fileID);


    end

end

end
