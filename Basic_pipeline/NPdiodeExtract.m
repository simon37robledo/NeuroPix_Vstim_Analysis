

% path = '\\sil3\data\\Large_scale_mapping_NP\lizards\PV139\PV139_Experiment_6_2_24\Insertion1\catgt_PV139_Experiment_6_2_24_1_g0';
%
% NP = NPAPRecording(path);

%%%%%Function inputs = NP, stimName, ismoving;

%stimName = "linearlyMovingBall";

%stimName = "rectGrid";

function [onsetSync offsetSync onSync offSync] = NPdiodeExtract(NP,ismoving)

%
%%%%MATLAB PTB SPECS: %%%%%%%%%%%%%%%%%%%%%%

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
end

%%%%%GET UN-SYNCED DIGITAL TRIGGERS %%%%%%%%%%%%%%%%%%%%%%


file = dir (NP.recordingDir);
filenames = {file.name};

Sstart = readmatrix(string(filenames(contains(filenames,"_tcat.nidq.xd_1_1_0"))))*1000;
Send = readmatrix(string(filenames(contains(filenames,"_tcat.nidq.xid_1_1_0"))))*1000;
onset = readmatrix(string(filenames(contains(filenames,"_tcat.nidq.xd_1_2_0"))))*1000;
offset = readmatrix(string(filenames(contains(filenames,"_tcat.nidq.xid_1_2_0"))))*1000;


ttl_index = 2;

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

FramesPsec = framesTrial/(stimDur/1000);

samplesPframe = 1/FramesPsec*(NP.samplingFrequencyNI);


for i = 1:length(stimOn)

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

    [pkvalsUp, pklocUp] = findpeaks(signal,'MinPeakProminence',2*stdS,'MinPeakDistance',samplesPframe*2-samplesPframe*0.1);

    [pkvalsDown, pklocDown] = findpeaks(fDat*-1,'MinPeakProminence',2*stdS,'MinPeakDistance',samplesPframe*2-samplesPframe*0.1);

    %filter out outlier peaks If there are more peaks than max

    if length(pklocUp) > nFrames/2-1 || length(pklocDown) > nFrames/2-1
        pklocUp = pklocUp(pkvalsUp<max(pkvalsUp)-std(pkvalsUp) & pkvalsUp>mean(pkvalsUp)-std(pkvalsUp));

        pklocDown = pklocDown(pkvalsDown<mean(pkvalsDown)+std(pkvalsDown) & pkvalsDown>mean(pkvalsDown)-std(pkvalsDown));
    end

    onFrame = [onFrame,pklocDown'./(NP.samplingFrequencyNI/1000)+stimOn(i)];

    onset = [onset min(pklocDown(1),pklocUp(1))./(NP.samplingFrequencyNI/1000)+stimOn(i)];

    offFrame = [offFrame, pklocUp'./(NP.samplingFrequencyNI/1000)+stimOn(i)];

    offset = [offset max(pklocUp(end),pklocDown(end))./(NP.samplingFrequencyNI/1000)+stimOn(i)];

end

% 
% figure()
% plot(AllD(startSnip:endSnip+500))
% xline((offset-stimOn(1))*NP.samplingFrequencyNI/1000)
% xline([pklocDown'])
% xline([pklocUp'],'red')% round(MBcrossDown(i)*)]+500)
% xline((stimOn(1)-MBcrossUp(i))*NP.samplingFrequencyNI/1000,'red')
% xline(500,'blue')
% xline((stimOff(1)-MBcrossUp(i))*NP.samplingFrequencyNI/1000,'red')
%
%
% xline(onFrame*(NP.samplingFrequencyNI/1000)+500,'green')
% xline(offFrame*(NP.samplingFrequencyNI/1000)+500,'red','.')
%
% figure()
% plot(AllD(round(length(AllD)-10*NP.samplingFrequencyNI):end))
% xline((onFrame(onFrame>-10000+MBcrossDown(end))-MBcrossUp(1))*NP.samplingFrequencyNI/1000-round(length(AllD)-10*NP.samplingFrequencyNI));
% xline((offFrame(offFrame>-10000+MBcrossDown(end))-MBcrossUp(1))*NP.samplingFrequencyNI/1000,'g')

else

    onset = [];
    offset = [];
    for i =1:length(stimOn)

        startSnip  = round((stimOn(i)-stimOn(1))*(NP.samplingFrequencyNI/1000))+1;
        endSnip  = round((stimOff(i)-stimOn(1)+500)*(NP.samplingFrequencyNI/1000));

        if endSnip>length(AllD)
            signal = AllD(startSnip:end);
            fDat=medfilt1(AllD(startSnip:end),2*frameSamples);
        else
            signal =(AllD(startSnip:endSnip));
            fDat=medfilt1(AllD(startSnip:endSnip-1),2*frameSamples);
        end

        stdS = std(fDat);

        [pkvalOn, pklocOn] = findpeaks(-1*fDat,'MinPeakProminence',2*stdS,'MinPeakDistance',(stimDur-100)*NP.samplingFrequencyNI/1000);
        
        pklocOn = min(pklocOn);
        afterOn = find(signal>mean(fDat));
        pklocOfflist = find(afterOn>pklocOn);
        pklocOff = afterOn(pklocOfflist(1));

        onset = [onset pklocOn./(NP.samplingFrequencyNI/1000)+stimOn(i)];
        offset = [offset pklocOff./(NP.samplingFrequencyNI/1000)+stimOn(i)];        
     

    end

    figure()
    plot(signal)
    hold on;plot(fDat)
    xline(pklocOn,'red')
    xline(pklocOff)
    yline(mean(fDat))
    xline((stimOn(1:5)-stimOn(1))*NP.samplingFrequencyNI/1000)

end

%%%% SYNC %%%%%%%%%%%%%%%%%%%%%%
fileName = "DIODE_"+string(stimName);


% Search in for neural SQW and NI SQW and load them
cd(NP.recordingDir)
Neur = readmatrix(dir(fullfile(NP.recordingDir, '*imec0.ap.xd_384_6_500.txt*')).name);
originSQW = readmatrix(dir(fullfile(NP.recordingDir, '*nidq.xd_1_6_500.txt*')).name);
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
offsetSync = interp1(originSQW'*1000, Neur'*1000, offset, 'linear');



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
