%%%Plot all receptive fields 


sign = 0.05;
for ex =  51%GoodRecordingsPV%SDGrecordingsA%GoodRecordings%GoodRecordingsPV%GoodRecordingsPV%selecN{1}(1,:) %1:size(data,1)
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

    RFuSTDir = load(sprintf('RFuSelecTimeD-%s',NP.recordingName)).RFuSTDir;

    RFu = load(sprintf('RFuSelecTime-%s',NP.recordingName)).RFuST;

    RFuNormVid = load(sprintf('NormVideo-%s',NP.recordingName)).NormVideo;

    p = NP.convertPhySorting2tIc(NP.recordingDir);
    label = string(p.label');
    goodU = p.ic(:,label == 'good');

    respNeuronsMB = load(sprintf('pvalsBaselineBoot-%d-%s',N_bootstrap,NP.recordingName)).pvalsResponse;

    respU = find(respNeuronsMB<sign);
    
    bin =1;

    [Mb] = BuildBurstMatrix(goodU(:,respU),round(p.t/bin),round((directimesSorted-preBase)/bin),round(preBase/bin));

    [nT,nN,nB] = size(Mb);
    Nb2=  mean(Mb,3);

    Nbase = mean(Mb,[1 3]);

    normMatrixMean =reshape(Nbase,[1,1,nN]);
    normMatrixSTD = reshape(std(Nb2),[1,1,nN])+eps;

    RFuNormVidSTD = sum(RFuNormVid,3).*normMatrixSTD;
    RFuNormVidMean = sum(RFuNormVid,3).*normMatrixMean;


    RFuNorm = (RFu - RFuNormVidMean)./RFuNormVidSTD;


  


    if eNeuron <= 10

        tiles1 = eNeuron;

        tiles2 = 1;

    else

        tiles1 = 10;
        tiles2 = ceil(length(eNeuron)/10);

    end

    

    eNeuron = 7;
    for u = eNeuron
        fig = tiledlayout(direcN/2,direcN/2);

        for d = 1:direcN
            if d ==1 || d==3
                nexttile;imagesc(rot90(squeeze(RFuSTDir(d,:,:,u)),2));c = colorbar;caxis([0 max(RFuSTDir(:,:,:,u),[],'all')]);
            else
                nexttile;imagesc((squeeze(RFuSTDir(d,:,:,u))));c = colorbar;caxis([0 max(RFuSTDir(:,:,:,u),[],'all')]);
            end
            colormap('jet')
            title(string(uDir(d)))
            title(c,'Z-score')
            %xlim([(redCoorX-redCoorY)/2 (redCoorX-redCoorY)/2+redCoorY])

%             xt = xticks;%(linspace((redCoorX-redCoorY)/2,(redCoorX-redCoorY)/2+redCoorY,offsetN*2));
%             xticklabels(round(theta_x(1,xt)))
%             yt = yticks;
%             yticklabels(round(theta_y(yt,1)))
%             xlabel('X degrees')
%             ylabel('Y degrees')

        end

    end

end