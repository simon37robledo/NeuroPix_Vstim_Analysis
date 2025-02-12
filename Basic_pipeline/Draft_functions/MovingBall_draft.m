

%% Find shape of receptive field and whether it is signifficanlty different than the rest of the image. 
fun = @(block_struct) ...
   mean(block_struct.data,'all');
for u = 1:size(goodU,2)
    for s = 1:sizeN

        M=squeeze(normRFu(s,:,:,u));
        M2 = blockproc(M,[5 5],fun);
        NaNs = isnan(M);
        M(isnan(M))=-1;
        M2 = blockproc(M,[5 5],fun);
        %M=M+1;

        figure;imagesc(M);colorbar
        M2n = M;%normalize(M2,'center','mean');

        maxBlock = max(M2n(~NaNs),[],'all');
        minBlock = min(M2n(~NaNs),[],'all');

        respDelta = sum(maxBlock-M2n((~NaNs)),'all',"omitnan")/numel(M2n((~NaNs)));
        RTI = 1-((respDelta)^2/(maxBlock-minBlock)^2);

    end

end


%%

% Load an image
img = squeeze(normRFu(:,:,17));

% Convert to grayscale if it's not
if size(img,3) == 3
    img = rgb2gray(img);
end

% Define block size
blockSize = [50, 50]; % example size

% Define the function to find maximum
fun = @(block_struct) max(block_struct.data(:));

% Apply the function to each block
maxBlock = blockproc(img, blockSize, fun);

% Display the result
imshow(maxBlock, []);

figure(normRFu(normRFu~=-1))
 td= squeeze(normRFu);
 td = td(td~=-1);
 figure;imagesc(td);

figure;imagesc(squeeze(normRFu(:,:,:,17)));colorbar
figure;imagesc(squeeze(RFu(:,:,:,17)));colorbar
figure;imagesc(squeeze(matrixResp(:,:,20)))
figure;imagesc(squeeze(matrixForNorm));colorbar
figure;imagesc(matrix_nXnY);

unique(matrix_nXnY)
unique(matrixResp)
unique(RFu(:,:,20))

colorbar


%% Calculate the spike sum per frame

spikeSum = zeros(nT,length(goodU),sizeX(4));
%spikeSumBase = mean(Mbase,3);

% Make raster matrix that contains only stimulus. 
msPerFarme= round(stimDur/sizeX(4));
[Mr] = BuildBurstMatrix(goodU,round(p.t),round((directimesSorted)),round((stimDur)));
[Mbase] = BuildBurstMatrix(goodU,round(p.t),round((directimesSorted-preBase)),round((preBase)));

for u = 1:length(goodU)

    Mu = squeeze(Mr(:,u,:));

    j= 1;  

    for f = 1:sizeX(4)  
        
        spikeSum(:,u,f) = sum(Mu(:,1*j:min(f*msPerFarme,length(Mu))),2);

        j = f*msPerFarme;
    end
end
spikeSum = single(spikeSum);
spikeRates=spikeSum./(msPerFarme/1000); %Divide by seconds per frame
    

%% Prepare for convolution
%%get videos of every trial
% Define the grid size
xAxis = 1:cell2mat(ball.VSMetaData.allPropVal(63));%63 for pv103,64 for pv67 %65 PV139 and on   %:max(X,[],'all')-min(X,[],'all')-max(sizeV)/2;

yAxis = 1:cell2mat(ball.VSMetaData.allPropVal(63));%63 PV103,64 for pv67 %65 and on   %max(Y,[],'all')-min(Y,[],'all')-max(sizeV)/2;

%respU = [9,11,12,13,15,16,17,18,19,20,21,22,23,24,25,31]; %PV139_2

respUgt =[18,21,22,24,26,30,31,32,34,35,38,39,41,42,44,47,50,51,53,57,58,61,63,65,66]; %PV139_1

%respU = [66, 65, 61, 42, 24,  21]; %PV59_1

%respU = [24]; %PV32_1
%respU = [2,8,16]; %PV32_2
%respU = [45,9]; %PV32_3
% respU = [18,22,24,26,42,46,47,50,53,56,61,64,65,68,73,79,81,84,89,90,94,97,98,100,101,102,103,104,107,108,112,...
%      113,119,120,121,124,127,130,131,133,137,138,141,142,144,145,147,149,151,154,157]; %PV67_1
%respU = [115,102,77,75]; %PV67_4
%respUgt = [1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 17, 20, 22, 23, 24, 28, 32, 34, 36]; %PV103_1
%respU = [2,23,25,28,20,22,25,28,29,32,34,36,37,38,39,40]; %PV103_2
%respU = [3,7,10,14,15,17,18,19,20,26,33,46,50]; %PV103_3
%respU = [5,8,13];%PV103_4
%respU = [35,89];%PV103_5
%respU = [11,13,14,21,29,31,32,28,45,80,81,82]; %PV103_6
%respU = [6,18,25,29,30,32,36,37,38,44,57,78,82,84,89]; %PV103_7
%respU = sort(respU);

cd(NP.recordingDir + "\Figs")

fixed_delay = 1;


RFuDir = cell(1,direcN);
RFuSize = cell(1,sizeN);
RFuSpeed = cell(1,speedN);

A = [stimOn directions' offsets' sizes' speeds'];

C = sortrows(A,[5 2 3 4]);

uDir = unique(directions);
uSize = unique(sizes);
uSpeed = unique(speeds);

for i = 1:direcN
    C(C(:,2) == uDir(i),2) = i;
end
for i = 1:sizeN
    C(C(:,4) == uSize(i),4) = i;
end
for i = 1:speedN
    C(C(:,5) == uSpeed(i),5) = i;
end

centerXreal = cell2mat(ball.VSMetaData.allPropVal(61)); %PV139 63 %PV64 62 %PV103 61 
%%%%%% Convolution
%respU =[144,145,147,149,151,154,157];
%%%Change conv so that the true geometry is used, and that it doesn't
%%%iterate through neurons. 

%respU = 18;
%for u = 1:length(respU)
    
    RFu = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');

    for k = 1:direcN
        RFuDir{k} = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');
    end
    for k = 1:sizeN
        RFuSize{k} = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');
    end
    for k = 1:speedN
       RFuSpeed{k} = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');
    end

    for i = 1:trialDivision:trials

        videoTrials = zeros(length(xAxis)+1,length(yAxis)+1,sizeX(4),'single');
        for j = 1:sizeX(4) %%Calculate video of unique trials

            % Create a matrix of zeros representing the grid
            xyScreen = zeros(length(xAxis)+1,length(yAxis)+1,'single');

            %Define the center and radius of the circle

%             centerX = ChangePosX(i,j)-min(X,[],'all'); % X position of the center
%             centerY = ChangePosY(i,j)-min(Y,[],'all'); % Y position of the center
            centerX = ChangePosX(i,j)-(centerXreal-length(yAxis)/2); % X position of the center
            centerY = ChangePosY(i,j); % Y position of the center

            radius = sizeV(i)/2; % Radius of the circle

            % Create a meshgrid of coordinates
            [x, y] = meshgrid(1:length(xAxis)+1, 1:length(yAxis)+1);

            % Calculate the distance of each point from the center
            distances = sqrt((x - centerX).^2 + (y - centerY).^2);

            % Set the values inside the circle to 1 (or any other value you prefer)
            xyScreen(distances <= radius) = 1;

            videoTrials(:,:,j) = xyScreen;
            %figure;imagesc(xyScreen);
        end

        spikeMean = mean(spikeSum(i:i+trialDivision-1,respU(u),:))-mean(Mbase(i:i+trialDivision-1,respU(u),:),'all'); %Get the mean across same trials and substract the baseline

        tic
        Co = convn(videoTrials,spikeMean,'same');
        toc
        %toc

%         spikeMean = mean(spikeSum(i:i+trialDivision-1,respU(u),:));
% 
%         tic
%         Co = convn(videoTrials,spikeMean,'same');
%         toc

        RFuDir{C(i,2)} =  RFuDir{C(i,2)}+Co;
        RFuSize{C(i,4)} = RFuSize{C(i,4)}+Co;
        RFuSpeed{C(i,5)} = RFuSpeed{C(i,5)}+Co;

        RFu = RFu+Co;
    end

    save(sprintf('RFu-%d',respU(u)),'RFu')
    save(sprintf('RFuDir-%d',respU(u)),'RFuDir','-v7.3')
    save(sprintf('RFuSize-%d',respU(u)),'RFuSize','-v7.3')
    save(sprintf('RFuSpeed-%d',respU(u)),'RFuSpeed','-v7.3')
%end
%%




%%
%Delay = 204 ms = 12 bins of 17 ms

delay = 12;
if rem(sizeX(4),2)==0
    center = sizeX(4)/2+1;
else
    center = (sizeX(4)-1)/2+1;
end

%%Search for maximal intensity
% figure
% plot([1:154],squeeze(spikeMean)')

%% Depth distribution of responses

basicPathA = {basic_pathPV67,basic_pathPV103,basic_pathPV27,basic_pathPV139,basic_pathPV59,basic_pathPV32};

expA = {expPV67,expPV103,expPVPV27,expPV139,expPV59,expPV32};

InsAnglesPV67 = [88 88 88 88];
InsDepthsPV67 = [3956, 3907, 4000, 3934]; 
RespUnitPV67 ={[18,22,24,26,42,46,47,50,53,56,61,64,65,68,73,79,81,84,89,90,94,97,98,100,101,102,103,104,107,108,112,...
    113,119,120,121,124,127,130,131,133,137,138,141,142,144,145,147,149,151,154,157],[],[],[115,102,77,75]};

InsAnglesPV103= [70.5 70.5 70.5 70.5 70.5 70.5 70.5 ];
InsDepthsPV103= [4104.14, 3964.18, 4066.5, 4123.87, 4175.86, 4225.55, 4027.42]; 
RespUnitPV103={[1, 2, 3, 7, 8, 9, 10, 11, 12, 13, 17, 20, 22, 23, 24, 28, 32, 34, 36],[2,23,25,28,20,22,25,28,29,32,34,36,37,38,39,40],...
    [3,7,10,14,15,17,18,19,20,26,33,46,50],[5,8,13],[35,89],[11,13,14,21,29,31,32,28,45,80,81,82],[6,18,25,29,30,32,36,37,38,44,57,78,82,84,89]};

InsAnglesPV27= [72.5 72.5 72.5  72.5 72.5 72.5 72.5];
InsDepthsPV27= [3913.9 3904.34 3525.7 3900.1 3914.34 3739.3 3906.8];
RespUnitPV27={[16,]};

InsAnglesPV139 = [89 89];
InsDepthsPV139 = [3910 3907.23];
RespUnitPV139 = {[18,21,22,24,26,30,31,32,34,35,38,39,41,42,44,47,50,51,53,57,58,61,63,65,66],...
    [9,11,12,13,15,16,17,18,19,20,21,22,23,24,25,31]};

InsAnglesPV59 = [81 81 81];
InsDepthsPV59 = [2845.89 3358.91 3050.36];
RespUnitPV59 = {[66, 65, 61, 42, 24,  21],[],[]};

InsAnglesPV32 = [69 72 72];
InsDepthsPV32 = [2400 2600 2200];
RespUnitPV32 = {[24],[2,8,16],[45]};

animalD = {InsDepthsPV67,InsDepthsPV103,InsDepthsPV27,InsDepthsPV139,InsDepthsPV59,InsDepthsPV32};

animalA = {InsAnglesPV67,InsAnglesPV103,InsAnglesPV27,InsAnglesPV139,InsAnglesPV59,InsAnglesPV32};

RespUnits = cell2mat({RespUnitPV103,RespUnitPV67,RespUnitPV139,RespUnitPV59,RespUnitPV32});

verticalDepth = cell(2,length(animalA));

YDist = cell(2,length(animalA));
%% Depth plot
figure()
j = 1;
aj=2;
for a = [1,2,4:length(animalA)] %X doesn't change


    %point1 = [animals{a}(:,1),animals{a}(:,2),animals{a}(:,3)];

    verticalDepth{1,a} = sin(deg2rad(animalA{a})).*(animalD{a}); %depth of unit along vertical axis

    YDist{1,a} = cos(deg2rad(animalA{a})).*(animalD{a}); %X distance of unit from insertion

    %point2 = [animals{a}(:,1),animals{a}(:,2)+YDist',animals{a}(:,3)-verticalDepth'];

    %     for i = 1:length(animals{a})
    %         plot3([point1(i,1),point2(i,1)],[point1(i,2),point2(i,2)],[point1(i,3),point2(i,3)],'LineWidth', 1, 'Color',colorA{a})
    %         hold on
    %     end

    for in = 1:length(animalA{a})
        path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
            in+string(filesep)+"catgt_"+string(expA{a})+"_"+in+"_g0");

        NP = NPAPRecording(path);

        cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t');

        %Good units

        GoodU_or = cluster_info.cluster_id(cluster_info.group=="good");

        GoodU_orDepth = cluster_info.depth(cluster_info.group=="good");

        %GoodU_orDepth = GoodU_orDepth(RespUnits{a}{in});

        verticalDepth{2,a}{in} = sin(deg2rad(animalA{a}(in)))*(animalD{a}(in) - GoodU_orDepth);

        ResponsiveDepth =  verticalDepth{2,a}{in}(RespUnits{a}{in});
        Nindex = true(1, length(verticalDepth{2,a}{in}));
        Nindex(RespUnits{a}{in}) = false;
        NonRespDepth =  verticalDepth{2,a}{in}(Nindex);

        hold on;
        plot(ones(length(ResponsiveDepth),1)+j+0.5-1+aj-2, ResponsiveDepth,'o','MarkerFaceColor','b')
        hold on; plot(ones(length(NonRespDepth),1)+j-1+aj-2, NonRespDepth,'o','MarkerFaceColor','r')

        j=j+3;
    end
aj = aj+5;
end

ylabel('Unit depth (um)')
xlabel('insertions and animals (blue= visually responive units to Moving ball stim)')


%% find unique pairs and their counts. Create Normalization matrix
%coords = [ChangePosX(:)-min(X,[],'all')+1,ChangePosY(:)-min(Y,[],'all')+1];
%centerXreal = cell2mat(ball.VSMetaData.allPropVal(61)); %63 in PV139 and on
%coords = [ChangePosX(:)+1-(centerXreal-length(yAxis)/2),ChangePosY(:)+1];%translate coordinates to 1080x1080

%task: ask Mark about exact parameters of screen, what is the appropiate parameter? 

difX = max(ChangePosX,[],'all') - min(ChangePosX,[],'all'); %max(X,[],'all')-min(X,[],'all')
difY = max(ChangePosY,[],'all') - min(ChangePosY,[],'all');
% 
% nX = max(coords(:,1));
% nY = max(coords(:,2));
%[Cn, ~, ic] = unique(coords, 'rows', 'stable');
%counts = accumarray(ic, 1);
find(arrayfun(@(x) strcmp(x.allPropName , 'rect'), ball.VSMetaData));
coorRect = cell2mat(ball.VSMetaData.allPropVal(57));
find(contains(ball.VSMetaData.allPropName,'rect'));%63 in PV139 and on

%[x, y] = meshgrid(1:length(xAxis)+1, 1:length(yAxis)+1);
[x, y] = meshgrid(1:coorRect(3),1:coorRect(4));
%[x, y] = meshgrid(1:difX,1:difY);
sizesU = unique(sizeV);

matrixForNorm = zeros(length(xAxis)+1, length(yAxis)+1);
%matrixForNorm  = zeros(1920+1,1920+1);
%matrixForNorm  = zeros(difX,difY);

matrixForNorm = zeros();

for s = 1:sizeN
    % Step 3: Fill in the matrix with the counts of repeated coordinates
    for i = 1:trialDivision*sizeN:length(ChangePosX)
        %matrix_nXnY = zeros(length(xAxis)+1, length(yAxis)+1);
        matrix_nXnY = zeros(coorRect(4),coorRect(3));
        %matrix_nXnY = zeros(difX,difY);
        for tr = 1:sizeX(4)

            %             centerX = Cn(i,1);
            %             centerY = Cn(i,2);

            centerX = ChangePosX(i,tr);%-(centerXreal-length(yAxis)/2);
            %centerX = ChangePosX(i,tr)*length(xAxis)/coorRect(3);
            centerY = ChangePosY(i,tr);
            %centerY = (ChangePosY(i,tr)+(centerXreal-length(yAxis)/2))*length(yAxis)/coorRect(3);

            % Create a meshgrid of coordinate
            radius = sizesU(s)/2;

            % Calculate the distance of each point from the center
            distances = sqrt((x - centerX).^2 + (y - centerY).^2);

            % Set the values inside the circle to 1 (or any other value you prefer)
            matrix_nXnY(distances <= radius) =1;
           
        end
        matrixForNorm = matrixForNorm+matrix_nXnY;
    end
      
end

figure;imagesc(matrixForNorm)
colorbar
%needs to be transformed into a 1080x180
%PV103 done / PV64 / PV139 
expA = {expPV67,expPV103,expPVPV27,expPV139,expPV59,expPV32};
cd(NP.recordingDir)
save(sprintf('matrixForNorm-%s',expA{4}),'matrixForNorm')


%%%Message for mark:
% Regarding the moving ball graphs, when I try to reconstruct a summary image of the stimulus, I get something weird. Either not all offsets are shown, or all offsets are there but some don't intersect, as if the extreme offsets where off the screen (I checked by running a stim with the same parameters, and indeed all offsets intersect). This is a problem of identifying on which coordinate reference the X and Y of the moving ball are recorded. I am using  the rect parameter which is a [0 0 1920 1080]. This results in offsets not intersecting. 
%   


%% Plot RF with a stationary delay
%to do:
%1. Standarize  the screen size to 1080 (for animal 6, bigger screen size)
%2. Detect maximum screen value per unit, and normalize, then sum across
%units. 
%3. Find better way to detect delay in convolution (max of video?). 
%4. Save result matrices. 


cd(NP.recordingDir+"\Figs")
%respU = 25;
% respU
% %respU = 157;
sideVF = 1080;%cell2mat(ball.VSMetaData.allPropVal(63)); %65 in pv139 and on

% SumAllDir = cell(2,length(animalA)); 
% SumPerDir = cell(2,length(animalA)); 
sizeScreenAn = [64,63,64,65,65,65]; %changes in the parameters order. 
%%
for a =[2]
    inSumAll = cell(1,length(animalA{a}));
    
    inSumDir = cell(2,length(animalA{a}));
    

    for in = 1:length(animalA{a})
        respU = RespUnits{a}{in};

        %generate path
        path = convertStringsToChars(string(basicPathA{a})+string(filesep)+string(expA{a})+string(filesep)+"Insertion"+...
            in+string(filesep)+"catgt_"+string(expA{a})+"_"+in+"_g0");
        NP = NPAPRecording(path);
        
        %Load stim stats
        patternIndex = strfind(string(NP.recordingDir), "\catgt");
        endIndex = patternIndex(1)-1;
        stimDir = string(NP.recordingDir);
        stimDir = extractBetween(stimDir,1,endIndex);
        file = dir (stimDir);
        filenames = {file.name};
        ballFiles = filenames(contains(filenames,"linearlyMovingBall"));
        ball= load(stimDir+"\"+string(ballFiles));
        directions = cell2mat(ball.VSMetaData.allPropVal(17));
        uDirections = unique(directions);
        sizeScreen = cell2mat(ball.VSMetaData.allPropVal(sizeScreenAn(a)));
        
        cd(path+"\Figs")

        Allu = zeros(sizeScreen+1,sizeScreen+1,length(respU));
        AlluDir = zeros(sizeScreen+1,sizeScreen+1,length(respU),length(uDirections));        

        for u =1:length(respU)
            
            %respU(u) =157;
            %sideVF = cell2mat(ball.VSMetaData.allPropVal(63)); %65 in pv139 and on

            M = load(sprintf('RFu-%d.mat',respU(u)));
            M = M.RFu;
            Mdir = load(sprintf('RFuDir-%d.mat',respU(u)));
            Mdir= Mdir.RFuDir;

            [sM1 sM2 sM3] = size(M);

            if sM1 >  sideVF+1 || sM2 >  sideVF+1

                leftPix = floor((sM1-sideVF)/2);

                rightPix = sM1-floor((sM1-sideVF)/2)-1;

                M = M(leftPix:rightPix,leftPix:rightPix,:);

            end

%             fig = figure; % Open a new figure window
%             imagesc(squeeze(M(:,:,center+delay)))%.
%             %axis equal; % Ensure the aspect ratio is equal to make the circle look r ound
%             %axis off; % Optionally, turn off the axis for a cleaner look
%             title(sprintf('U.%d-ReceptiveField-Del-%d ms',respU(u),delay*msPerFarme));
%             set(gca,'YDir','normal')
%             colorbar
%             %prettify_plot
%             cd(path+"\Figs")
%             %print(fig, sprintf('%s-MovBall-U%d-RFconv.png',NP.recordingName,respU(u)),'-dpng');
%             close

            Allu(:,:,u) = M(:,:,center+delay);

          %  t = tiledlayout(length(Mdir)/2,length(Mdir)/2);

            for direc = 1:length(Mdir)

                %nexttile; % Open a new figure window
                Md = Mdir{direc};
               % imagesc(squeeze(Md(:,:,center+delay)))%.
                %axis equal; % Ensure the aspect ratio is equal to make the circle look r ound
                %axis off; % Optionally, turn off the axis for a cleaner look
                %title(sprintf('U.%d-ReceptiveField-Del-%d ms',respU(u),delay*msPerFarme));
                %set(gca,'YDir','normal')
                %colorbar
                %caxis([0, 0.6]);
                %prettify_plot

                AlluDir(:,:,u,direc) =  Md(:,:,center+delay);
                %dirsSum{2,direc} = uDirections(direc);
                
            end     
%             cd(path+"\Figs")
%             print(fig, sprintf('%s-MovBall-U%d-RFconv.png',NP.recordingName,respU(u)),'-dpng');
%             close

        end
        inSumAll{in} =  Allu;
        inSumDir{1,in} = AlluDir;
        inSumDir{2,in} = uDirections;
    end
    SumAllDir{1,a} = inSumAll;
    SumAllDir{2,a} = expA{a};
    SumPerDir{1,a}= inSumDir;
    SumPerDir{2,a}= expA{a};
end


%% Get summary images. 


for a = [1,2,4]

    f1= figure;
    subplot(1,3,1)

    for in=1:length(animalA{a})

        uA = SumAllDir{1,a}{in};
        meanA = mean(uA,3);

        %imagesc(meanA-


    end

end
%%



M = load('RFu-26.mat');
M = M.RFu;
[aaa, aaaain] = max(M,[],3);

[aa,aain] = max(aaa,[],'all');

aaaain(aain); 
fig = figure; % Open a new figure window
imagesc(squeeze(M(:,:,center+delay)))%./matrixForNorm*spkRateBM(respU));
%colormap('gray'); % Use a grayscale colormap





M = load('RFu26.mat');
M = M.RFu;

% %% Test conv center
% spikeMean(1,1,:) = zeros(1,1,154); 
% spikeMean(1,1,78) = 1; %--> center = 78 154/2+1
% C = convn(videoTrials,spikeMean,'same');



% %%Video File Test
% %M = rand(240, 320, 100);
% videoFile = 'u29.avi';
% writerObj = VideoWriter(videoFile);
% writerObj.FrameRate = 60; % Set the frame rate to 30 frames per second
% open(writerObj);
% 
% for k = 1:size(M, 3) % Loop through each frame
%     frame = M(:, :, k); % Extract the k-th frame
%     % If your data isn't uint8, you might need to normalize and convert it
%     %frame = im2uint8(mat2gray(frame)); % Normalize and convert to uint8
%     writeVideo(writerObj, frame); % Write the frame to the video
% end
% 
%    close(writerObj);
% 
%    implay(videoFile);

%%
    
% %%Iterate per neuron
% 
% for u = 1:length(goodU)
% 
%     Mu = squeeze(Mr(:,u,:));
% 
%     j= 1;
% 
%     spikeSum = zeros(trials,sizeX(4));
% 
%     for f = 1:sizeX(4)  
%         
%         spikeSum(:,f) = sum(Mu(:,1*j:min(f*msPerFarme,length(Mu))),2);
% 
%         j = f*msPerFarme;
% 
%     end
% 
%     iter = max(spikeSum,[],'all'); %iteration number
% 
%     for i = 1:iter
% 
%         [row,col] = find(spikeSum>=i);
% 
%         Ysum = ChangePosY(row,col);
%         Xsum = ChangePosX(row,col);
% 
%         X_flat = Xsum(:); % Flatten A into a column vector
%         Y_flat = Ysum(:); % Flatten B into a column vector
% 
%         % Combine into pairs
%         pairs = [X_flat, Y_flat];
% 
%         % Find unique pairs and their counts
%         [unique_pairs, ~, ic] = unique(pairs, 'rows', 'stable');
%         pair_counts = accumarray(ic, 1);
% 
%         %Transoform pair intro screen coordinates
%         x = round(unique_pairs(:,1))-round(min(X,[],'all'))+1;
%         y = round(unique_pairs(:,2))-round(min(Y,[],'all'))+1;
% 
%         for c = 1:length(pair_counts)
%             xyScreen(y(c), x(c)) = pair_counts(c);
%         end
% 
%         nGrid = 5;
%         
%         redGrid = zeros(nGrid,nGrid);
% 
%         for b = 1:nGrid
% 
%             for v = 1:nGrid
% 
%                 div = round(length(xyScreen)/nGrid);
% 
%                 redGrid(b,v) = mean(xyScreen(b*div-div+1:b*div,v*div-div+1:v*div),'all');
% 
%             end
% 
%         end
%         figure
%         
%         imagesc(redGrid);
% 
%         max(y)
% 
% 
%     end
% 
% 
%     for x = 1:length(xAxis)
%         for y = 1:lrngth(yAxis)
% 
%             xpos=find
% 
%         end
%     end
% 
% 
% 
% end
% 


