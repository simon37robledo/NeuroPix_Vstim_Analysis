

% path = '\\sil3\data\\Large_scale_mapping_NP\lizards\PV139\PV139_Experiment_6_2_24\Insertion1\catgt_PV139_Experiment_6_2_24_1_g0';
%
% NP = NPAPRecording(path);

%%%%%Function inputs = NP, stimName, ismoving;

%stimName = "MovBall";

%stimName = "rectGrid";

function [onsetSync offsetSync onSync offSync] = NPdiodeExtract(NP,ismoving,stimName,ttl_index,digCH,syncCH) %dig -16/1 sync 0/6

%
%%%%MATLAB PTB SPECS: %%%%%%%%%%%%%%%%%%%%%%
cd(NP.recordingDir)

fileName = "DIODE_"+string(stimName);


if exist(NP.recordingDir+"\"+fileName+"OnFrame"+".txt",'file') & exist(NP.recordingDir+"\"+fileName+"OffFrame"+".txt",'file')...
        & exist(NP.recordingDir+"\"+fileName+"Onset"+".txt",'file') & exist(NP.recordingDir+"\"+fileName+"Offset"+".txt",'file') & ismoving

    onsetSync = readmatrix(fileName+"Onset"+".txt");
    offsetSync =  readmatrix(fileName+"Offset"+".txt");
    onSync = readmatrix(fileName+"OnFrame"+".txt");
    offSync = readmatrix(fileName+"OffFrame"+".txt");

elseif exist(NP.recordingDir+"\"+fileName+"Onset"+".txt",'file') & exist(NP.recordingDir+"\"+fileName+"Offset"+".txt",'file') 

    onsetSync = readmatrix(fileName+"Onset"+".txt");
    offsetSync =  readmatrix(fileName+"Offset"+".txt");

else


    patternIndex = strfind(string(NP.recordingDir), "\catgt");

    endIndex = patternIndex(1)-1;
    stimDir = string(NP.recordingDir);
    stimDir = extractBetween(stimDir,1,endIndex);

    file = dir (stimDir);
    filenames = {file.name};

    if ismoving


        ballFiles = filenames(contains(filenames,"linearlyMovingBall"));

        directions = [];
        offsets = [];

        j =1;
        if size(ballFiles) ~= [0 0]

            for i = ballFiles
                ball= load(stimDir+"\"+string(i));

                j = j+1;
            end
            disp('Visual stats extracted!')
        else
            disp('Directory does not exist!');
        end
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

    stimDur = stimOff(1) - stimOn(1);%transform into horizontal vector.


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

        nFrames = unique(cell2mat(ball.VSMetaData.allPropVal(21)));

        framesTrial = ball.VSMetaData.allPropVal{21,1}(1,1,1);

        FramesPsec = 60;%ball.VSMetaData.allPropVal{60,1};

        samplesPframe = 1/FramesPsec*(NP.samplingFrequencyNI);


        for i = 1:length(stimOn)


            pklocUpN = [];
            pklocDownN = [];
            %i=44;
            startSnip  = round((stimOn(i)-stimOn(1))*(NP.samplingFrequencyNI/1000))+1;
            endSnip  = round((stimOff(i)-stimOn(1)+100)*(NP.samplingFrequencyNI/1000));

            if endSnip>length(AllD)
                signal = AllD(startSnip:end);
                fDat=medfilt1(AllD(startSnip:end),15);
            else
                signal =(AllD(startSnip:endSnip));
                fDat=medfilt1(AllD(startSnip:endSnip-1),15);
            end

            stdS = std(signal);

            [pkvalsUp, pklocUp] = findpeaks(signal,'MinPeakProminence',0.8*stdS,'MinPeakDistance',samplesPframe*2-samplesPframe*0.3);

            [pkvalsDown, pklocDown] = findpeaks(fDat*-1,'MinPeakProminence',0.8*stdS,'MinPeakDistance',samplesPframe*2-samplesPframe*0.3);

            %filter out outlier peaks If there are more peaks than max
            %

            pklocUp = pklocUp(pklocUp>1100);
            dF = diff(pklocUp); %Filter out initial false peaks contained in a trial that starts late. 
            dFf = find(dF<mean(dF)*2);
            pklocUp = pklocUp(dFf(1):end);

            pklocDown = pklocDown(pklocDown>1100); %% 100 is the normal delay between the difital trigger and the analog onset. 

%             if length(pklocUp) > nFrames/2 || length(pklocDown) > nFrames/2
%                 pklocUp = pklocUp(pkvalsUp<max(pkvalsUp)-0.5*std(pkvalsUp)&kvalsUp>max(pkvalsUp)+0.5*std(pkvalsUp));
% 
%                 pklocDown = pklocDown(pkvalsDown<mean(pkvalsDown)+0.5*std(pkvalsDown)&pkvalsDown>mean(pkvalsDown)-0.5*std(pkvalsDown));
%             end

            %             firstF = min(pklocDown(1),pklocUp(1));
            %             if mean(signal(1:firstF))>mean(signal)
            %                 iFrame = find(signal(1:firstF)<mean(signal(1:firstF)));
            %             else
            %                 iFrame = find(signal(1:firstF)>mean(signal(1:firstF)));
            %             end
            %
            %             iniFrame = iFrame(1);
            %
            %             lastF = max(pklocDown(end),pklocUp(end));
            %             if mean(signal(lastF:end))>mean(signal)
            %                 eFrame = find(signal(lastF:end)<mean(signal(lastF:end)));
            %             else
            %                 eFrame = find(signal(lastF:end)>median(signal(lastF:end)));
            %             end
            %
            %             endFrame = lastF+eFrame(end);
            %
            %             eframe = pklocUp+eFrame(1);

            %%%Start frame:

            if (length(pklocUp) < nFrames/2 || length(pklocDown) < nFrames/2) && pklocUp(1) < pklocDown(1)

                pklocDownN = [pklocUp(1)-samplesPframe;pklocDown];

            elseif (length(pklocUp) < nFrames/2 || length(pklocDown) < nFrames/2) && pklocUp(1) > pklocDown(1)

                pklocUpN = [pklocDown(1)-samplesPframe;pklocUp];

            end
            %%%%end frame:
            if (length(pklocUp) < nFrames/2 || length(pklocDown) < nFrames/2) && pklocUp(end) < pklocDown(end)

                pklocUpN = [pklocUp;pklocDown(end)+samplesPframe;];

            elseif (length(pklocUp) < nFrames/2 || length(pklocDown) < nFrames/2) && pklocUp(end) > pklocDown(end)

                pklocDownN = [pklocDown;pklocUp(end)+samplesPframe];

            end

            if isempty(pklocUpN)
                pklocUpN = pklocUp;
            end
            if isempty(pklocDownN)
                pklocDownN = pklocDown;
            end

            %%%%%%%

%             if isempty(pklocDownN)||isempty(pklocUpN)
% 2+2
%             end

            [startP indS] =  min([pklocDownN(1),pklocUpN(1)]);

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



        end
% 
%         %
%         figure()
%         plot(signal)
%         xline([pklocDown'])
%         xline([pklocUp'],'red')% round(MBcrossDown(i)*)]+500)
%         %         xline(pklocDown(1)-16.66*NP.samplingFrequencyNI/1000,'blue')
        %         xline(endFrame,'green')

        %         xline((stimOn(1)-MBcrossUp(i))*NP.samplingFrequencyNI/1000,'red')
        %         xline(500,'blue')
        %         xline((stimOff(1)-MBcrossUp(i))*NP.samplingFrequencyNI/1000,'red')
        %
        %
        %         xline(onFrame*(NP.samplingFrequencyNI/1000)+500,'green')
        %         xline(offFrame*(NP.samplingFrequencyNI/1000)+500,'red','.')
        %
        %         figure()
        %         plot(AllD(round(length(AllD)-10*NP.samplingFrequencyNI):end))
        %         xline((onFrame(onFrame>-10000+MBcrossDown(end))-MBcrossUp(1))*NP.samplingFrequencyNI/1000-round(length(AllD)-10*NP.samplingFrequencyNI));
        %         xline((offFrame(offFrame>-10000+MBcrossDown(end))-MBcrossUp(1))*NP.samplingFrequencyNI/1000,'g')

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

end

end
