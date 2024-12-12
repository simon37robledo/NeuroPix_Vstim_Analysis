%%% Test moving window RG

bin=50;
win=stimDur+stimInter;
preBase = round(stimInter/20)*10;

[Mr]=BuildBurstMatrix(goodU,round(p.t/bin),round((directimesSorted-preBase)/bin),round(win/bin));

[nT,nN,nB] = size(Mr);
%indRG --> sorted infexes

trialsPerCath = length(seqMatrix)/(length(unique(seqMatrix)));

PositionsTotal = positionsMatrix(seqMatrix(indRG),:);


[posS,indexX] = sortrows(PositionsTotal,1); %Sort first dimension because tile layout moves through columns

MrS = Mr(indexX,:,:);

posX = squeeze(NeuronVals(:,:,3));
posY = squeeze(NeuronVals(:,:,2));

for u =1;%1:length(goodU)
    t = tiledlayout(sqrt(max(seqMatrix)), sqrt(max(seqMatrix)),'TileSpacing','tight');



    if trialsPerCath>20
        mergeTrials = 5;

        Mr2 = zeros(nT/mergeTrials,nB);

        for i = 1:mergeTrials:nT

            meanb = mean(squeeze(MrS(i:min(i+mergeTrials-1, end),u,:)),1);

            Mr2(j,:) = meanb;

            j = j+1;

        end
    else
        Mr2=MrS(:,u,:);
        mergeTrials =1;
    end


    [T,B] = size(Mr2);
    j=1;
    for i = 1:trialsPerCath/mergeTrials:T
        %Build raster
        M = Mr2(i:min(i+trialsPerCath/mergeTrials-1, end),:);
        [nTrials,nTimes]=size(M);
        nexttile
        imagesc((1:nTimes),1:nTrials,squeeze(M));colormap(flipud(gray(64)));
        xline(preBase/bin, '--', LineWidth=2, Color="#77AC30");
        xline((stimDur+preBase)/bin, '--', LineWidth=2, Color="#0072BD");
        xticks([preBase/bin (round(stimDur/100)*100+preBase)/bin]);
        if nSize >1
            yline([trialsPerCath/mergeTrials/nSize:trialsPerCath/mergeTrials/nSize:trialsPerCath/mergeTrials-1]+0.5,LineWidth=1)
        end
        caxis([0 1]);
        set(gca,'YTickLabel',[]);

        if i < T - (trialsPerCath/mergeTrials)*max(positionsMatrix(:))-1
            set(gca,'XTickLabel',[]);

        end
    
        hold on
        rectangle('Position', [round(posX(u,j)/bin)+round(preBase/bin), 0, round(duration/bin), trialsPerCath/mergeTrials+1],...
            'EdgeColor', 'r', 'LineWidth', 1.5,'LineStyle','-.');

        j = j+1;
    end
    fig = gcf;
    set(fig, 'Color', 'w');
    colorbar

    %Plot rectangle:

 
    
    % Set the color of the figure and axes to black
    title(t,sprintf('Rect-GRid-raster-U%d',u))
    ylabel(t,sprintf('%d trials',nTrials*mergeTrials))
    fig.Position = [227         191        1413         781];
    %prettify_plot
    print(fig,sprintf('%s-rect-GRid-raster-U%d.png',NP.recordingName,u),'-dpng')
    close
end


