cd('\\sil3\data\Large_scale_mapping_NP')
excelFile = 'Experiment_Excel.xlsx';

data = readtable(excelFile);

%%
for ex = 40%examplesSDG%[7 8 28]%1:size(data,1)
    %%%%%%%%%%%% Load data and data paremeters
    %1. Load NP class
    path = convertStringsToChars(string(data.Base_path(ex))+filesep+string(data.Exp_name(ex))+filesep+"Insertion"+string(data.Insertion(ex))...
        +filesep+"catgt_"+string(data.Exp_name(ex))+"_"+string(data.Insertion(ex))+"_g0");
    try %%In case it is not run in Vstim computer, which has drives mapped differently
        cd(path)
    catch
        originP = cell2mat(extractBetween(path,"\\","\Large_scale"));
        if strcmp(originP,'sil3\data')
            path = replaceBetween(path,"","\Large_scale","W:");
        else
            path = replaceBetween(path,"","\Large_scale","Y:");
        end
        cd(path)
    end
    NP = NPAPRecording(path);

    %Create Figs and matData folders if they don't exist

    if ~exist(path+"\Figs",'dir')
        mkdir Figs
    end

    if ~exist(path+"\matData",'dir')
        mkdir matData
    end
    cd('C:\Users\MarkS9\Documents\GitHub')

 %2. Extract moving ball statistics
    patternIndex = strfind(string(NP.recordingDir), "\catgt");

    endIndex = patternIndex(1)-1;
    stimDir = string(NP.recordingDir);
    stimDir = extractBetween(stimDir,1,endIndex);

    file = dir (stimDir);
    filenames = {file.name};

    file = dir (stimDir);
    filenames = {file.name};
    ballFiles = filenames(contains(filenames,"linearlyMovingBouncingDots"));

    
        if isempty(ballFiles)
            %disp()
            w= sprintf('No moving ball files where found in %s. Skipping into next experiment.',NP.recordingName);
            warning(w)
            continue
        end


    dotNumbers = [];
    dotSize = [];
    statictime = [];
    speeds = [];
    xi = [];
    vi = [];

    j =1;


    if size(ballFiles) ~= [0 0]

        for i = ballFiles %Extract stim properties
            ball= load(stimDir+"\"+string(i));

            
            dotNumbers = [dotNumbers cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'dotsNumbers'))))];

            dotSize = [dotSize cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'dotSize'))))];
            
            statictime = [statictime cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'waitFromOnset'))))];

            speeds = [speeds cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'allSpeeds'))))];

            xi = [xi ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'xi')))];

            vi = [vi ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'vi')))];

            interStimStats = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            stimDurSt = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'stimDuration'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

 
    %static_images

    imgSize = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));
    imgSize = sort(imgSize(3:4));

%%
    grayLevel = 0.5;  % Gray background (range 0 to 1, where 0=black, 1=white)
    % Initialize the image with gray

    % Create a meshgrid for the image coordinates
    %[cols, rows] = meshgrid(1:imgSize(2), 1:imgSize(1));
    [cols, rows] = meshgrid(1:round(imgSize(2)/dotSize), 1:round(imgSize(1)/dotSize));
    cols = single(cols);
    rows = single(rows);
    nTrials = length(xi{1});

    % Draw each dot on the image

    Images = cell(1,nTrials);

    for tr = 1:nTrials
        x = single(round(xi{1}{tr}(:,1)./dotSize)); 
        y = single(round(xi{1}{tr}(:,2)./dotSize));
        image = zeros(round(imgSize./dotSize),'single');
        for i = 1:dotNumbers
            % Compute distance from the current dot's center
            distance = sqrt((cols - x(i)).^2 + (rows - y(i)).^2);
            
            % Set pixels within the radius of the dot to white
            image(distance <= 1/2) = 1;  % White dot
        end
        Images{tr} = image;

    end

    figure;imagesc(image);
    figure;imagesc(Images{1});
    %%
%
    ifi = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ifi'))));
    nFrames = ceil((stimDurSt/1000)/ifi);
    visualFieldRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'visualFieldRect'))));
    visualFieldRect = round(visualFieldRect./dotSize);
    videosTrial = cell(nTrials,1);
    cd(NP.recordingDir)
    %save('videosTrialBouncingB.mat',"videosTrial",'-v7.3')
    %[cols, rows] = meshgrid(1:round(imgSize(2)/dotSize), 1:round(imgSize(1)/dotSize));
    %videosTrial = matfile('videosTrialBouncingB.mat','Writable',true);

    if ~isfile('videosTrialBouncingB.mat')
        for i =1:nTrials
            x=single(xi{1}{i});
            v=single(round(vi{1}{i}./dotSize));
            xAll=zeros([size(x),nFrames],'single');
            xAll(:,:,1)=round(x./dotSize);
            vidTrial = zeros(size(cols,1),size(cols,2),nFrames,'single');
            for f=2:nFrames
                %grayLevel = 0.5;  % Gray background (range 0 to 1, where 0=black, 1=white)
                image = zeros(round(imgSize./dotSize),'single');
                xAll(:,:,f)=xAll(:,:,f-1)+v*ifi;
                p1=xAll(:,1,f)>=visualFieldRect(3) | xAll(:,1,f)<=visualFieldRect(1);
                p2=xAll(:,2,f)>=visualFieldRect(4) | xAll(:,2,f)<=visualFieldRect(2);
                pC=p1|p2;
                xAll(pC,:,f)=xAll(pC,:,f-1)-v(pC,:)*ifi;
                v(p1,1)=-v(p1,1);
                v(p2,2)=-v(p2,2);
                xAll(pC,:,f)=xAll(pC,:,f-1)+v(pC,:)*ifi;

                xCoor = xAll(:,1,f);
                yCoor = xAll(:,2,f);

                for d = 1:dotNumbers
                    % Compute distance from the current dot's center
                    distance = sqrt((cols - xCoor(d)).^2 + (rows - yCoor(d)).^2);
                    % Set pixels within the radius of the dot to white
                    image(distance <= 1/2) = 1;  % White dot
                end
                vidTrial(:,:,f) = image;
            end
            vidTrial(:,:,1) = round(Images{i});
            videosTrial{i} = vidTrial;
            %videosTrial.videosTrial(i,1) = vidTrial;
        end
        save('videosTrialBouncingB.mat',"videosTrial",'-v7.3')
    else
        videosTrial = load('videosTrialBouncingB.mat').videosTrial;
    end
    %implay(vidTrial)
    


    %%
    %3. Load Triggers (diode)
    Ordered_stims= strsplit(data.VS_ordered{ex},',');
    containsBB = cellfun(@(x) contains(x,'BB'),Ordered_stims);
    ttlInd = find(containsBB);

    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,0,1,"BB",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex));
    [stimOn stimOff onSync offSync] = NPdiodeExtract(NP,0,1,"BB",ttlInd,data.Digital_channel(ex),data.Sync_bit(ex)); %Ugly second time to make sure orientation is right for creating A


    stimInter= mean(stimOn(2:end)-stimOff(1:end-1)); % When dealing with different speeds, save different stim durs, and create M for each speed
    stimDur = mean(-stimOn+stimOff); %When dealing with different speeds, save different stim durs, and create M for each speed
    
    
    %5. Load data to select good units
    cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');
    GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");
    GoodU_orDepth = cluster_info.depth(cluster_info.group=="good");

    %If new tic matrix needs to be used, move oldTIc files and run convertPhySorting2tIc
    cd(NP.recordingDir)

    p = NP.convertPhySorting2tIc(NP.recordingDir);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');

    staticTime = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'waitFromOnset'))))*1000;

    delay = 200;
    bin = ifi*1000;
    [Mr] = BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn'+staticTime+delay)/bin),nFrames); %response matrix of onset of movement
    [nT,nN,nB] = size(Mr);

    RFu = zeros(nN,size(videosTrial{1},1),size(videosTrial{1},2));
%
    for u =1:length(goodU)%1:nN
        for tr = 1:nT

            Mrtr = (Mr(tr,u,:)).*videosTrial{tr};
            RFu(u,:,:) = squeeze(RFu(u,:,:))+sum(Mrtr,3);

        end
    end

figure;imagesc(squeeze(RFu(u,:,:)));

%%% Occupation matrix
OccM = zeros(size(videosTrial{1},1),size(videosTrial{1},2));

for tr = 1:nT
    OccM = OccM + sum(videosTrial{tr},3);
end



figure;imagesc(OccM);

%%% Normalized by occupation:

%baseline:
[Mb] = BuildBurstMatrix(goodU,round(p.t/bin),round((stimOn'-interStimStats)/bin),round((nFrames*ifi*1000+2000)/bin));

MbUnit = squeeze(mean(squeeze(mean(Mb)),2));

%Norm matrix per neurons:
NormUnit = zeros(nN,size(videosTrial{1},1),size(videosTrial{1},2));

for u = 1:nN
    NormUnit(u,:,:) = OccM.*MbUnit(u);
end


%%%% To convolve with a circular kernel
% Create a circular kernel (disk) with a radius of 5
radius = 5;
se = strel('disk', radius, 0);  % Structuring element (disk shape)
circularKernel = se.Neighborhood;  % Convert the structuring element to a binary matrix

% Normalize the circular kernel by dividing it by its sum
normalizedKernel = circularKernel / sum(circularKernel(:));
meanResult = zeros(nN,size(videosTrial{1},1),size(videosTrial{1},2));

% Apply 2D convolution with 'same' size option to keep the original matrix size
for u = 1:length(goodU)
    meanResult(u,:,:) = conv2(squeeze(RFu(u,:,:)), normalizedKernel, 'same');
    NormUnit(u,:,:) = conv2(squeeze(NormUnit(u,:,:)), normalizedKernel, 'same');
end

figure;imagesc(squeeze(meanResult(u,:,:)));
title(sprintf('delay=%d',delay))

figure;imagesc(squeeze(NormUnit(u,:,:)));
title(sprintf('delay=%d',delay))

%%% Normalized:
NormR = meanResult./NormUnit;

figure;imagesc(squeeze(NormR(31,:,visualFieldRect(1):visualFieldRect(3))));
title(sprintf('delay=%d',delay))



%%

end