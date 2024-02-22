%% Receptive field moving ball

% TO DO: understand code
%0. Check baseline calculation and Z/score, use non stim as baseline. 
%1. Implement response delay parameter (around 200 ms). 
%1. Use diode times if there is a good diode signal, other wise use
%estimation. 
%2. Make any directions available.

%inputs -> NP, stimDir, heatmap, Nrand units, Nresponsive units
%%% 
%Call example:

dgTR = Ttrigger{5};

aaM = ReceptiveField_MovStim(NP,p,stimDir, stimOn, stimDur,stimInter, RespU, 1);

%%%

function [ Mrf ] = ReceptiveField_MovStim(NP, p, stimDir, stimOn, stimDur, stimInter, respU, varargin)

%%% Function calculates receptive fields for moving stimuli and can plot it in
%%% a heat map (1 per neuron). The output is a cell array containing an x
%%% by x matrix M per neuron, where x is the number of cells in the grid
%%% representing the screen. 

%%MANDATORY INUTS:
%- NP --> neuropixels class
%- stimDir --> directory of the stimulus parameters from Psychotoolbox
%- respU --> Responsive units accordying to previous statistic. 
%%OPTIONAL INPUTS:
%- heatmap --> 0 no heatmpa is plotted | 1 heatmap is plotted
%- NrandR --> Number of responsive units sampled randomly.If 0 then plot
%all responsive neurons
%- NrandNR --> Number of non-responsive units sampled randomly. If 0 then plot
%all non-responsive neurons
%- specifU --> Vector of specific units that are to be plotted.

%%%%Optional inputs defaults>
heatmap = 0;
NrandR = 1;
NrandNR = 1;
specifU = 0;
lenGrid =7;
preBase = 500;
DiodeGood = 1;

%nargin
% Update optional parameters if provided


if nargin > 7
    heatmap = varargin{1};
end
if nargin > 8
    NrandR = varargin{2};
end
if nargin > 9
    NrandNR = varargin{3};
end
if nargin > 10
    specifU = varargin{4};
end
if nargin > 11
    lenGrid = varargin{5};
end
if nargin >12
    preBase = varargin{6};
end
if nargin >12
    DiodeGood = varargin{7};
end


cluster_info = readtable(string(NP.recordingDir) + "\cluster_info.tsv",  "FileType","text",'Delimiter', '\t'); %Get KS output

%Good units

GoodU_or = cluster_info.cluster_id(cluster_info.group=="good"); %Select good units

if specifU ~= 0
    GoodU_or = specifU; %Check if units are specified in input to function
end

%Select good units from p matrix 
label = string(p.label');
goodU = p.ic(:,label == 'good');

%Initiate output array 
Mrf = zeros(length(GoodU_or), lenGrid, lenGrid);

% Get Matlab PTB data

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

%To order offsets and directions:
A = [stimOn' directions' offsets'];

C = sortrows(A,[2 3]);

directimesSorted = C(:,1)';

%Check if there are more than 4 directions:

if direcNames ~= [0 1.5707963267949 3.14159265358979 4.71238898038469]

    directions = directions(ismember(directions,direcNames));

end

%Load spike index and spike times:

spkI = readNPY(string(NP.recordingDir) + "\spike_clusters.npy");

spkt = readNPY(string(NP.recordingDir) + "\spike_times.npy")/(NP.samplingFrequency/1000);

%%%

%Load X and Y positions from PTB (organized by offsets, directions, and
%ball positions:

X = squeeze(cell2mat(ball.VSMetaData.allPropVal(22)));

Y = squeeze(cell2mat(ball.VSMetaData.allPropVal(23)));

%Get the times values of frames per trial based on the stimulus duration and
%the length of X (option without diode):

times = linspace(0,round(stimDur),length(X));

sX = size(X);

%Get number of trials per offset-direction cathegory:

trials = cell2mat(ball.VSMetaData.allPropVal(29));

%Change order of dimensions of positions X & Y, so the reshaped vector
% order is frames-Ofsset1-Dir1, then frames-offset2-dir1.... then
% frames-offset1-dir2:

Xp =  permute(X, [3 1 2]); 

Xv =  reshape(Xp, 1, []); %%Position vector.

n = numel(Xv);

%Get number of positions per direction:
direcDiv = length(Xv)/length(unique(directions));

%
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


[nT,nN,nB] = size(M);

%%%2.Convolute in the 1st dimension (bin) and reduce
%%%element number to number of offsets

trialDivision = nT/(length(unique(offsets))*length(unique(directions)));


stdBL = std(M(:,:,1:round(baseLine/frameDur)),0,3);

meanBL = mean(M(:,:,1:round(baseLine/frameDur)),3);

stdBLperN = std(meanBL);

meanBLperN = mean(meanBL);


%replace zero values

stdBL(stdBL==0) = mean(stdBLperN);

stdBLperN(stdBLperN==0) = mean(stdBLperN);

allM1 = zeros(length(GoodU_or),lenGrid,lenGrid);
allM2 = zeros(length(GoodU_or),lenGrid,lenGrid);
allM3 = zeros(length(GoodU_or),lenGrid,lenGrid);
allM4 = zeros(length(GoodU_or),lenGrid,lenGrid);


%%%ITERATE ACROS UNITS

for u = 1:length(GoodU_or)

    u = 116;

    %change channel 0 to 1

    cluster_info.ch(cluster_info.ch ==0) = 1;
    ch = cluster_info.ch(cluster_info.cluster_id == GoodU_or(u));
    title_stim = sprintf('Unit-%d-channel-#%d', GoodU_or(u), ch);

    %%%Build raster:
    bin = 30;
    win = stimDur+preBase*2;
    [Mraster]=BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-baseLine)/bin),round(win/bin));

    %Now to reduce the number of trials so spikes can be seen:
    redYDim = 4;
    filter = (1/redYDim)*ones(4,1);
    MRred =conv2(squeeze(Mraster(:,u,:)),filter,'valid');
    MRred = MRred(redYDim:redYDim:end,:);
    MRred(MRred>0) = 1;

    [nT,nB] = size(MRred);
    offsetN = length(unique(offsets));
    direcN = length(unique(directions));
    TrialsPerOffset = nT/(offsetN*direcN);

    fig = figure;
    imagesc(MRred);colormap(flipud(gray(64))); %plot it
    xline(preBase/bin,'k', LineWidth=1.5) %Stim start
    xline(stimDur/bin+preBase/bin,'k',LineWidth=1.5) %stim end
    ylabel('Trials');xlabel('Time (ms)');
    title(title_stim);
    xticks([0.5 (preBase/bin):10:nB])
    xticklabels([-preBase 0:10*bin:nB*bin])
    DirectionLines = nT/direcN:nT/direcN:nT-1;
    yline(DirectionLines+0.5,'r--', LineWidth=3);
    yline(TrialsPerOffset:TrialsPerOffset:nT-1,'k',LineWidth=0.5,Alpha = 0.2)
    yticklabels(yticks*redYDim)
    caxis([0 1])
    prettify_plot

    %%%%

    %Main unit channel position

    Unit_spks = spkt(spkI == GoodU_or(u));


    sC =[]; %%%Spike times

    for i= 1:length(t1) %iterate every time and asign number of spikes

        sC = [sC sum(Unit_spks > t1(i) & Unit_spks < t1(i)+(t1(2)-t1(1)))/(frameDur)];
    end

    BL = zeros(1,length(t1));

    for i=1:sX(3):length(t1) %Iterate to get baseline (1000 ms before start)

        BL(i:i+sX(3)-1) = sum((Unit_spks > t1(i)-interStimStats/4 & Unit_spks < t1(i)))/(baseLine);

        %sC(i:i+sX(3)-1) = (sC(i:i+sX(3)-1) - BL(i:i+sX(3)-1))/std(BL(i:i+sX(3)));

    end

    B = [t1;rx;ry;sC;BL]; %Create matrix

    %%%hERE ITERATE, until here the code is general

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

    Rms = 100; %response period in ms

    frameMs = stimDur/length(times);

    framResp = round(Rms/frameMs);

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
    nM1 = M1;%((M1-mean(concatMb,'all',"omitnan"))/stdBLrect);

    nM2 = M2;%((M2-mean(concatMb,'all',"omitnan"))/stdBLrect);

    nM3 = M3;%((M3-mean(concatMb,'all',"omitnan"))/stdBLrect);

    nM4 = M4;%((M4-mean(concatMb,'all',"omitnan"))/stdBLrect);


    nM1 = nM1(2:end-1,2:end-1);
    nM2 = nM2(2:end-1,2:end-1);
    nM3 = nM3(2:end-1,2:end-1);
    nM4 = nM4(2:end-1,2:end-1);

    %%%GENERAL
    nM = ((M-mean(concatMb,'all',"omitnan"))/stdBLrect);
    %%%%

    maxVal = unique(max(cat(2, nM1,nM2,nM3,nM4),[],'all'));
    minVal =  unique(min(cat(2, nM1,nM2,nM3,nM4),[],'all'));
    concatMn = cat(2, nM1,nM2,nM3,nM4);

    if heatmap == 1

        %%Plot

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
        %caxis([-max(abs(concatMn(:)))-0.1 max(abs(concatMn(:)))]);
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
        %set(gca,'XTick',[],'YTick',[]);
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
        %set(gca,'XTick',[],'YTick',[]);
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
        %set(gca,'XTick',[],'YTick',[]);
        %ax.CLim = [0  maxVal];
        %colorbar('gray');

        title(t,sprintf('%s-MovBall-Positions.png',title_stim))

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


