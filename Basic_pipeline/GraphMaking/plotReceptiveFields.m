%%%Plot all receptive fields 


% fig = tiledlayout(length(GoodRecordingsPV),1);
sign = 0.05;
j = 1;

%43 PV35_4 is empty

for ex = GoodRecordingsRF%SDGrecordingsA%GoodRecordings%GoodRecordingsPV%GoodRecordingsPV%selecN{1}(1,:) %1:size(data,1)
    
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


    %%%%%%%%%%Load stim statistics:


    %2. Extract moving ball statistics
    patternIndex = strfind(string(NP.recordingDir), "\catgt");

    endIndex = patternIndex(1)-1;
    stimDir = string(NP.recordingDir);
    stimDir = extractBetween(stimDir,1,endIndex);

    file = dir (stimDir);
    filenames = {file.name};

    file = dir (stimDir);
    filenames = {file.name};
    ballFiles = filenames(contains(filenames,"linearlyMovingBall"));

    
        if isempty(ballFiles)
            %disp()
            w= sprintf('No moving ball files where found in %s. Skipping into next experiment.',NP.recordingName);
            warning(w)
            continue
        end

    directions = [];
    offsets = [];
    sizes = [];
    speeds = [];
    orientations = [];

    j =1;

    if size(ballFiles) ~= [0 0]

        for i = ballFiles %Extract stim properties
            ball= load(stimDir+"\"+string(i));

            if ~isempty(find(strcmp(ball.VSMetaData.allPropName,'orientations'))) %Check if orientations are present (grid inside moving object).
                orientations = [orientations cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'orientations'))))];
            else %No orientations property exist
                orientations = [orientations zeros(1,cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials')))))];
            end
            if isempty(orientations) %orientations property exists but is empty
                orientations = [orientations zeros(1,cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'nTotTrials')))))];
            end

            directions = [directions cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'directions'))))];

            offsets = [offsets cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'offsets'))))];

            sizes = [sizes cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'ballSizes'))))];

            speeds = [speeds cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'speeds'))))];

            interStimStats = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'interTrialDelay'))))*1000;

            j = j+1;
        end
        disp('Visual stats extracted!')
    else
        disp('Directory does not exist!');
    end

    if ~isempty(find(strcmp(ball.VSMetaData.allPropName,'orientationFreq')))
        Freq = [ cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'orientationFreq'))))];
    end

    uDir = unique(directions);

    direcN = length(uDir);

    coorRect = cell2mat(ball.VSMetaData.allPropVal(find(strcmp(ball.VSMetaData.allPropName,'rect'))));

    sizeV = sizes(speeds == max(speeds))';
    reduceFactor = min([20 min(sizeV)]); %has to be bigger than the smallest ball size

    redCoorX = round(coorRect(3)/reduceFactor);
    redCoorY = round(coorRect(4)/reduceFactor);
   
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%Load Receptive fields %%%%%%%%%%%%%%%%

    RFuSTDir = load(sprintf('RFuSelecTimeD-%s',NP.recordingName)).RFuSTDir;

    RFuSTDirFilt = load(sprintf('RFuSTDirFilt-%s',NP.recordingName)).RFuSTDirFilt;

    RFuSTDirSizeFilt = load(sprintf('RFuSTDirSizeFilt-%s',NP.recordingName)).RFuSTDirSizeFilt; %Size and dir

    if size(RFuSTDirSizeFilt,2) >1 %If there is more than on size

        %%Select size closest to 120
        [minS indx] = min(abs(unique(sizeV)-120));

        RFuSTDirFilt = squeeze(RFuSTDirSize(:,indx,:,:,:));

    else

        RFuSTDirFilt = squeeze(RFuSTDirSize);

    end

%     RFu = load(sprintf('RFuSelecTime-%s',NP.recordingName)).RFuST;
% 
%     RFuNormVid = load(sprintf('NormVideo-%s',NP.recordingName)).NormVideo;
% 
%     p = NP.convertPhySorting2tIc(NP.recordingDir);
%     label = string(p.label');
%     goodU = p.ic(:,label == 'good');
% 
    respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;

    respU = find(respNeuronsMB<sign);

    %respU = 40;
    
%     bin =1;
% 
%     [Mb] = BuildBurstMatrix(goodU(:,respU),round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin));
% 
%     [nT,nN,nB] = size(Mb);
%     Nb2=  mean(Mb,3);
% 
%     Nbase = mean(Mb,[1 3]);
% 
%     normMatrixMean =reshape(Nbase,[1,1,nN]);
%     normMatrixSTD = reshape(std(Nb2),[1,1,nN])+eps;
% 
%     RFuNormVidSTD = sum(RFuNormVid,3).*normMatrixSTD;
%     RFuNormVidMean = sum(RFuNormVid,3).*normMatrixMean;
% 
%     RFuNorm = (RFu - RFuNormVidMean)./RFuNormVidSTD;


    if length(respU) <= 10

        tiles1 = length(respU);

        tiles2 = 1;

    else

        tiles1 = 10;
        tiles2 = ceil(length(respU)/10);

    end

    if numel(respU) ==1
        tiles1 = 1;
        tiles2 =1;


    end


%
    % %%%%%%Plot receptive fields
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fig = figure('Units', 'normalized', 'OuterPosition', [0 0 1 1]); % Full screen figure;
    InsertionLayout = tiledlayout(tiles1,tiles2,"TileSpacing","compact","Padding","none");
    %InsertionLayout.Layout.Tile = j;

    
    
    for u = 1:length(respU) %%%%problem with filtering. 

      
        %DirLayout = tiledlayout(direcN/2,direcN/2,"TileSpacing","tight");
        DirLayout=tiledlayout(InsertionLayout,direcN/2,direcN/2,"TileSpacing","tight","Padding","none");
        DirLayout.Layout.Tile = u;
       
        for d = 1:direcN
            
           nexttile(DirLayout);
            
            if d ==1 || d==3
               imagesc(rot90(squeeze(RFuSTDirFilt(d,:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY,u)),2));caxis([0 max(RFuSTDirFilt(:,:,:,u),[],'all')]);%c = colorbar;
            else
                imagesc((squeeze(RFuSTDirFilt(d,:,1+(redCoorX-redCoorY)/2:(redCoorX-redCoorY)/2+redCoorY,u))));caxis([0 max(RFuSTDirFilt(:,:,:,u),[],'all')]);%c = colorbar;
            end
            colormap('turbo')
            %title(string(uDir(d)))
            %title(c,'Z-score')
            %xlim([(redCoorX-redCoorY)/2 (redCoorX-redCoorY)/2+redCoorY])
            %xt = xticks;%(linspace((redCoorX-redCoorY)/2,(redCoorX-redCoorY)/2+redCoorY,offsetN*2));
            %xticklabels(round(theta_x(1,xt)))
            %yt = yticks;
            %yticklabels(round(theta_y(yt,1)))
            %             xlabel('X degrees')
            %             ylabel('Y degrees')

            xticks([])
            yticks([])

            axis equal tight
            

        end
       
        title(DirLayout,string(respU(u)),'FontSize',6)
       

    end

    title(InsertionLayout,NP.recordingName)
    cd('\\sil3\data\Large_scale_mapping_NP\Figs paper\Receptive fields all neurons\')
   

   % fig.Position = [ 0.506770833333333         0.037962962962963                 0.0734375         0.883333333333333];
    

    print(gcf, sprintf('Filt_Dir-RFs-%s.pdf',strrep(NP.recordingName,'_','-')), '-dpdf', '-r300', '-vector','-fillpage');

    j = j+1;

    close


end